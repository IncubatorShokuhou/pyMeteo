# 发布与贡献

[English](Publishing)

## 用户安装

```bash
pip install pymeteo-kit
```

导入仍为 `import pymeteo`。Python 3.6+。运行时只有 NumPy。

## PyPI

发行名：[`pymeteo-kit`](https://pypi.org/project/pymeteo-kit/)。发布走 GitHub Actions：

* 触发：GitHub **Release**，或 **Actions → Publish → Run workflow**
* 工作流：`.github/workflows/publish.yml`
* 环境：`pypi`（URL https://pypi.org/p/pymeteo-kit）
* 密钥：`PYPI_API_TOKEN`
* 构建：Python 3.12，`python -m build`（hatchling ≥ 1.18）

wheel/sdist 把 `core-metadata-version` 钉在 `2.1`，好让 CPython 3.6 上的 pip 21.3 仍能读元数据。

## CI

`.github/workflows/ci.yml`：3.12 上跑 ruff；3.8 / 3.10 / 3.12 / 3.13 上跑 pytest；3.6 job 把 3.12 构建的 wheel 装进 `python:3.6.15-buster`。

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

## 贡献（简）

* 问题与 PR 对着 https://github.com/IncubatorShokuhou/pyMeteo 的 `master`
* 包根公开 API 保持英文蛇形名。NCL 名字只放在 `pymeteo.ncl`；NCL 没有的名字不要伪造（热指数、风寒、沙氏、K、SWEAT 等）
* 单位用字符串。不要加 Pint、MetPy 或额外运行时依赖
* 只做点上诊断：除非经过明确评审的扩围，否则不做 I/O、绘图、地图投影、FAO56 全套、网格平流或完整 CAPE/CIN 探空套件
* 公式按公开文献自行实现，不要粘贴 MetPy 或 NCL 源码
* 改到的函数要有测试。NCL 文档例题数值放在 `tests/test_ncl_official_examples.py`

## 许可

MIT。版权人 IncubatorShokuhou，2019–2026。见源码树中的 `LICENSE`。
