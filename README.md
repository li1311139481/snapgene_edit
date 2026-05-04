
<div align="center">

# snapgene_edit

**纯 Python — 零外部依赖 — 读取、编辑、写入 SnapGene .dna 文件**

[![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)](#)
[![License](https://img.shields.io/badge/license-MIT-green)](#)

</div>

---

## 安装

### 方式一：直接从 GitHub 安装（推荐）

```bash
pip install git+https://github.com/li1311139481/snapgene_edit.git
```

### 方式二：clone 后以可编辑模式安装

```bash
git clone https://github.com/li1311139481/snapgene_edit.git
cd snapgene_edit
pip install -e .
```

> 可编辑模式修改源码后即时生效，适合开发和调试。

### 方式三：下载 ZIP 本地安装

```bash
# 先下载 ZIP 并解压
pip install path/to/snapgene_edit
```

---

## 快速验证

安装后运行以下代码确认成功：

```python
from snapgene_edit import SnapGeneEditor, batch_from_lists_simple

print("snapgene_edit 安装成功!")
```

---

## 基本用法

### 1. 读取 .dna 文件

```python
from snapgene_edit import SnapGeneEditor

editor = SnapGeneEditor("plasmid.dna")
print(f"名称: {editor.name}")
print(f"长度: {editor.length} bp")
print(f"环状: {editor.is_circular}")
print(f"特征数: {len(editor.features)}")

for f in editor.list_features():
    print(f"  - {f['name']} ({f['type']}) [{f['start']}:{f['end']}]")
```

### 2. 批量编辑 sgRNA 特征（克隆实验最常用）

```python
from snapgene_edit import batch_from_lists_simple

batch_from_lists_simple(
    save_names=["ZFP001", "ZFP002", "ZFP003", "ZFP004", "ZFP005", "ZFP006"],
    new_names=["sgKlf2-1", "sgKlf2-2", "sgKlf2-3", "sgKlf2-4", "sgKlf2-5", "sgKlf2-6"],
    new_seqs=[
        "GGGGGGGGGGGGGGGGGGGG",
        "AAAAAAAAAAAAAAAAAAAA",
        "TATATATATATATATATATA",
        "GCGCGCGCGCGCGCGCGCGC",
        "ACACACACACACACACACAC",
        "GGGGGGGGGGGGGGGGGGGG",
    ],
    fluor_list=["A", "B", "C", "G", "B", "B"],
    output_dir=r"/cluster/facility/hlhuang/zifeng001/1file/snapgene/",
    fluor_files={
        "A": r"/cluster/facility/hlhuang/zifeng001/1file/snapgene/LMA模板.dna",
        "B": r"/cluster/facility/hlhuang/zifeng001/1file/snapgene/LMB模板.dna",
        "C": r"/cluster/facility/hlhuang/zifeng001/1file/snapgene/LMC模板.dna",
        "G": r"/cluster/facility/hlhuang/zifeng001/1file/snapgene/LMG模板.dna",
    },
)
```

> `"FLUOR"` 会自动替换为对应荧光（A/B/C/G）的原始 sgRNA 序列。
> 输出包含处理日志和引物序列（正向 + 反向互补），可直接粘贴到 primer list。

---

## 配置

荧光模板文件路径在 `snapgene_edit/config.py` 中配置：

```python
FLUORESCENT_FILES = {
    "A": "D:/plasmid/LMA模板.dna",
    "B": "D:/plasmid/LMB模板.dna",
    "C": "D:/plasmid/LMC模板.dna",
    "G": "D:/plasmid/LMG模板.dna",
}
```

每个模板 .dna 文件中必须包含一个名为 `sgRNA` 的特征。

也可通过 `batch_from_lists_simple(..., fluor_files={...})` 参数传入自定义路径，无需修改配置文件。

---

## 项目结构

```
snapgene_edit/
├── pyproject.toml
├── example.py
├── snapgene_edit/
│   ├── __init__.py       # 公开 API
│   ├── editor.py         # 核心编辑器 + 批量函数
│   ├── reader.py         # .dna 读取
│   ├── writer.py         # .dna 写入
│   ├── models.py         # 数据模型
│   ├── config.py         # 荧光模板路径配置
│   └── utils.py          # 序列工具
└── docs/
    └── index.html        # 完整 API 文档
```

---

## 完整文档

详细 API 文档请查看 [docs/index.html](docs/index.html)（本地打开），或参考 [docs/index.html](https://github.com/li1311139481/snapgene_edit/blob/main/docs/index.html)。
