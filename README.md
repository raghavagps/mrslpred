# **MRSLpred**
A computational tool for multilabel mRNA subcellular localization prediction
## Introduction
MRSLpred is a tool for multilabel mRNA subcellular localization prediction using a XGBoost classifier. It uses only composition based features for predicting the subcellular locations. The final model also deploys a motif-based module which has been implemented using MERCI.
MRSLpred is also available as web-server at https://webs.iiitd.edu.in/raghava/mrslpred. Please read/cite the content about MRSLpred for complete information including algorithm behind the approach.

## Standalone
The Standalone version of mrslpred is written in python3 and following libraries are necessary for the successful run:
- scikit-learn
- xgboost=0.90
- Pandas
- Numpy


## Minimum USAGE
To know about the available option for the standlone, type the following command:
```
python3 mrslpred_motif.py -h
```
To run the example, type the following command:
```
python3 mrslpred_motif.py -i example_input.fa
```
This will predict where the submitted sequences are going to localize. It will use other parameters by default. It will save the output in "final_prediction.csv" and "final_prob_prediction.csv" in CSV format (comma separated variables).

## Full Usage
```
usage: mrslpred_motif.py [-h] --file FILE [--th1 TH1] [--th2 TH2] [--th3 TH3]
                         [--th4 TH4] [--th5 TH5] [--th6 TH6] --output OUTPUT

```
```
Please provide following arguments for successful run

optional arguments:
  -h, --help            show this help message and exit
  --file FILE, -f FILE  Path to fasta file
  --th1 TH1, -t1 TH1    Threshold for Ribosome
  --th2 TH2, -t2 TH2    Threshold for Cytosol
  --th3 TH3, -t3 TH3    Threshold for Endoplasmic Reticulum
  --th4 TH4, -t4 TH4    Threshold for Membrane
  --th5 TH5, -t5 TH5    Threshold for Nucleus
  --th6 TH6, -t6 TH6    Threshold for Exosome
  --output OUTPUT, -o OUTPUT    Path to output

```

**Input File:** It allow users to provide input in FASTA format.

**Output File:** Program will save the results to this folder

**Threshold 1:** User should provide threshold for Ribosome between 0 and 1, by default it is 0.3079.

**Threshold 2:** User should provide threshold for Cytosol between 0 and 1, by default it is 0.1468.

**Threshold 3:** User should provide threshold for Endoplasmic Reticulum between 0 and 1, by default it is 0.1156.

**Threshold 4:** User should provide threshold for Membrane between 0 and 1, by default it is 0.1956.

**Threshold 5:** User should provide threshold for Nucleus between 0 and 1, by default it is 0.7028.

**Threshold 6:** User should provide threshold for Exosome between 0 and 1, by default it is 0.9961.

MRSLpred Package Files
=======================
It contantain following files, brief description of these files given below

INSTALLATION                    : Installations instructions

LICENSE                         : License information

README.md                       : This file provide information about this package

xgboost_final.pkl               : This file contains the pickled version of model

mrslpred_motif.py               : Main python program

example_input.fa                : Example file contain nucleotide sequences in FASTA format

example_predict_prob_output.csv : Example output file containing probabilities for each location

example_predict_output.csv      : Example output file containing labels for each location
