# pChem

A modification-centric assessment tool for performance of chemoproteomic probes

# Introduction 

Chemical probe coupled with mass spectrometry (MS)-based proteomics, herein termed chemoproteomics, offers versatile tools to globally profile protein features and to systematically interrogate the mode of action of small molecules in a native biological system. Nonetheless, development of an efficient and selective probe for chemoproteomics can still be challenging. Besides, it is also difficult to unbiasedly assess its chemoselectivity at a proteome-wide scale. Here we present pChem, a modification-centric blind search and summarization tool to provide a pipeline for rapid and unbiased assessing of the performance of ABPP and metabolic labeling probes. This pipeline starts experimentally by isotopic coding of PDMs, which can be automatically recognized, paired, and accurately reported by pChem, further allowing users to score the profiling efficiency, modification-homogeneity and proteome-wide residue selectivity of a chemoproteomic probe.

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


We can then run the CMD as: 

```
python close_identify.py 
```

For more detailed usage, please refer to the instruction. 


