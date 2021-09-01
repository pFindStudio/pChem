# pChem

A modification-centric assessment tool for performance of chemoproteomic probes

# Introduction

This repository is an official implementation for **pChem**, which can generate all unknown modification candidates using blind search based on the pFind platform. 

For the ease of usage, we provide executable files (.exe) that can run in any environment without any additional package installation. Meantime, to improve the reproducibility and foster new researches in this field, we publicly release the source code of pChem. 

Please send any questions, comments or bug reports to feizhengcong@ict.ac.cn.  

# Requirements 

For exe version: 

* new version of pFind 

For python version: 

* new version of pFind 
* Python 3 
* matplotlib 
* pandas 

# Fast Usage 

In this section, we provide a simple and fast usage instruction to familiarize with basic operations. 

1. Set the parameter file *pChem.cfg*. You should set at least the path of fasta and msms dataset. Note that for pChem, we support the RAW or MZML formation for msms data. 

```
# ------------------------------------------------------------------
# parameter setting
# path to the pFind
pfind_install_path=D:\pFind3
# path to the fasta dataset
fasta_path=...
# path to the msms dataset
msmspath = ...
# ------------------------------------------------------------------
```

2. run the blind search module to find all the unknown modification candidates. 

In this module, we also make a accurate mass calibration, amino acid selectivity, and neural/charged losses recognition. 

We can run the CMD as: 

```
python blind_search.py 
```

3. run the close search module to make a restrcted search to further improve the resolution of dataset and identification efficiency of target PDM. 

```
python close_identify.py 
```

For more detailed usage, please refer to the instruction. 














