# 安装

[English](Install)

## 运行时

Python **3.6+**。运行时依赖只有 NumPy。CPython 3.6 上的 NumPy 上限是 **1.19.x**（最后一条完整支持 3.6 的发行线）。

```bash
pip install pymeteo-kit
```

PyPI 发行名是 `pymeteo-kit`（`pymeteo` 已被占用，`py-meteo` 因过于相似被拒）。导入仍为：

```python
import pymeteo as pm
```

## 从 git 检出

当前 `[build-system]` 要求 `hatchling>=1.18`，因此从源码做可编辑安装或 `python -m build` 需要 **Python 3.8+**。不要用 3.6 / 3.7 构建 sdist；这两个解释器请安装 PyPI 上的 `py3-none-any` wheel。

```bash
pip install -e .
pip install -e ".[dev]"   # pytest；ruff 在 3.8+
pytest
ruff check src tests
```

CI 里 lint 只跑 Python 3.12。本仓库 ruff 的 `target-version` 最低是 `py37`。

## 关于 3.6 的 CI

GitHub 托管的 `ubuntu-20.04` 与 `actions/setup-python` 的 3.6 镜像已于 2025 年下线。CI 仍用 `python:3.6.15-buster` 容器、对 3.12 构建的 wheel 跑 3.6 测试。wheel/sdist 把 `core-metadata-version` 钉在 `2.1`，好让 3.6 上最后一版 pip（21.3）能读元数据。

## 发布

GitHub Release，或手动 **Actions → Publish → Run workflow**，会构建 sdist/wheel 并用 `pypi` 环境里的 `PYPI_API_TOKEN` 上传到 PyPI。工作流文件：`.github/workflows/publish.yml`。详见 [发布与贡献](Publishing_zh)。
