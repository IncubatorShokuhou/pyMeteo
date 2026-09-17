# 热力学

[English](Thermo.md)

水汽和气块温度。库里大部分热力学计算在这一页。

函数可以从包根导入（`import pymeteo as pm`）。没特别说明时，默认单位是 °C、hPa、%、kg/kg。各函数下面写了用的公式。

## saturation_vapor_pressure

水面饱和水汽压。不单独处理冰面。

`temperature` 默认 `C`，结果默认 `hPa`。0 °C 时是 6.1078 hPa。

```python
import pymeteo as pm
pm.saturation_vapor_pressure(0.0)
pm.saturation_vapor_pressure(0.0, output_pressure_unit="Pa")  # 610.78
```

## condensation_temperature

抬升凝结高度上的温度，李社宏 1994 迭代。

参数是起始 `pressure`、`temperature`、`dewpoint`。网格插值若给出 Td > T 或负混合比，迭代可能停不下来，这时返回最后一次结果。

```python
pm.condensation_temperature(850.0, 16.6, 0.6)
```

## relative_humidity_from_dewpoint

温度和露点 → 相对湿度。输出默认 `%`；要 0–1 就设 `output_humidity_unit="fraction"`。

## dewpoint_from_relative_humidity

反过来。相对湿度 ≤ 0 的位置是 `nan`。

```python
pm.dewpoint_from_relative_humidity(18.0, 46.5)  # 大约 6.30 °C
```

## relative_humidity_from_mixing_ratio

温度、混合比、气压 → RH。饱和水汽压在 173.16–375.16 K 查表线性内插（和 NCL `relhum` 同一张表）。允许大于 100%；小于 0 截成 0.0001%。

`mixing_ratio_unit` 默认 `kg/kg`，也认 `g/kg`。

## mixing_ratio_from_relative_humidity

气压、温度、RH → 混合比。Tetens，对应 NCL `mixhum_ptrh`。

```python
pm.mixing_ratio_from_relative_humidity(1000.0, 18.0, 46.5)  # 约 0.00602 kg/kg
```

## specific_humidity_from_relative_humidity

参数一样，然后 `q = w / (1 + w)`。1000 hPa、18 °C、46.5% 时，`output_humidity_unit="g/kg"` 大约是 5.98 g/kg。

## convert_humidity

混合比和比湿互转，同时可以换质量单位。

`from_quantity` / `to_quantity` 取 `"mixing_ratio"` 或 `"specific_humidity"`。比湿 ≥ 1 再转混合比会得到 `nan`。

## visibility

用 RH 和温度估能见度。`method="RUC"` 是 60 km 乘指数衰减；`"FSL"` 用露点差。默认输出 km。

80%、18 °C、RUC 大约 11.8 km。

## saturation_mixing_ratio

`w_s = ε e_s(T) / (p − e_s)`，ε = 0.622。`e_s` 仍是李社宏水面公式。`e_s ≥ p` 时该点为 `nan`。

## mixing_ratio_from_dewpoint

把温度换成露点，公式相同。对应 NCL `mixhum_ptd` 那条路。

## vapor_pressure_from_mixing_ratio

`e = w p / (ε + w)`。混合比定义的代数逆。

## vapor_pressure_from_relative_humidity

`e = RH · e_s(T)`，`e_s` 还是李社宏公式。

## potential_temperature

Poisson 位温：`θ = T (1000 hPa / p)^0.286`。输出默认是 **K**，不是 °C。

```python
pm.potential_temperature(1000.0, 28.1)  # 301.25 K
```

## equivalent_potential_temperature

Bolton（1980）式 (43)，里面含他式 (22) 的 LCL 温度。输出默认 K。

## virtual_temperature

`T_v = T (1 + r/ε) / (1 + r)`。输出单位默认跟 `temperature_unit` 走。这是湿空气状态方程的精确形式，不用 NCL 文档里的 `T (1 + 0.61 r)` 近似。

## wet_bulb_temperature

Stull（2011）海平面经验湿球。大约按 1013 hPa 想的。比较靠谱的范围大致是 −20–50 °C、相对湿度 5–99%；又干又冷时偏差会大。

```python
pm.wet_bulb_temperature(20.0, 50.0)  # 大约 13.7 °C
```

## lifting_condensation_level

返回 `(p_lcl, T_lcl)`。`T_L` 用 Bolton（1980）式 (22)，再 `p_L = p (T_L / T)^{1/κ}`，κ = 0.286。

Wallace & Hobbs 那个 1000 hPa、15 °C、露点 4 °C 的例子，这里算出来大约 847 hPa。NCL `lclvl` 用 Stipanuk（1973），差几个 hPa 是正常的。

## parcel_temperature_at_pressure

气块先干绝热抬到 LCL，再湿绝热抬到 `pressure_target`。若目标气压还高于 LCL（`p_target ≥ p_LCL`），全程干绝热。`lifted_index_from_surface` 用的就是它。
