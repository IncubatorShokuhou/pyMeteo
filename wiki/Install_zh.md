# 安装

[English](Install)

Python **3.7+**。运行时依赖只有 NumPy。

```bash
pip install pymeteo-kit
```

PyPI 发行名是 `pymeteo-kit`（`pymeteo` 已被占用）。导入：

```python
import pymeteo as pm
```

从 git 检出：

```bash
pip install -e ".[dev]"
pytest
```

发布走 GitHub Actions 上传 PyPI。见 [发布与贡献](Publishing_zh)。
