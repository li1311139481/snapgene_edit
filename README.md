# snapgene_edit
读取、编辑和写入 SnapGene .dna 文件的纯 Python 库，无需 SnapGene 软件。目前仅支持序列长度相等的替换
# Tutorial
```
from snapgene_edit import batch_from_lists_simple

batch_from_lists_simple(
    save_names=["ZFP001", "ZFP002", "ZFP003", "ZFP004", "ZFP005", "ZFP006"],
    new_names=["sg1", "sg2", "sg3", "sg4", "sg5", "sg6"],
    new_seqs=[
        "GGGGGGGGGGGGGGGGGGGG",
        "AAAAAAAAAAAAAAAAAAAA",
        "TATATATATATATATATATA",
        "GCGCGCGCGCGCGCGCGCGC",
        "ACACACACACACACACACAC",
        "GGGGGGGGGGGGGGGGGGGG",
    ],
    fluor_list=["A", "B", "C", "G", "B", "B"],
    output_dir=r"D:/1lab_document/plasmid/克隆载体/",
)
```
