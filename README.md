# RAFBL

## Contents
- [RAFBL](#rafbl)
  - [Contents](#contents)
  - [About ↑](#about-)
  - [Install ↑](#install-)
  - [Run ↑](#run-)
    - [Generate features from Gaussian output files:](#generate-features-from-gaussian-output-files)
    - [Run whole model selection pipeline:](#run-whole-model-selection-pipeline)
    - [Screen for candidate ligands for the OA reaction:](#screen-for-candidate-ligands-for-the-oa-reaction)

## About [↑](#about) 
RAFBL is the repository accompanying the manuscript: _Reaction-Agnostic Featurization of Bidentate Ligands for Bayesian Ridge Regression of Enantioselectivity_. It includes two packages `modsel` and `moltop`.

`modsel` is used for additional ligand featurization from base features and takes care of the feature selection for the final models.

`moltop` generates topological features from molecular structures. A molecular graph is either constructed using xyz coordinates and covalent radii or SMILES directly.  

## Install [↑](#install) 
We recommend the use of [conda](https://docs.conda.io/en/latest/miniconda.html) to install all the require dependencies.

To create the environment, run:
`conda env create -f environment.yml`
And then activate the environment as:
`conda activate rafbl`

## Run [↑](#run) 

### Generate features from Gaussian output files:

To re-generate the features from Gaussian log files you can run:
```bash
./feat_csd.sh
./feat_lit.sh
```
This process takes a long time but only has to be run once. Beware! If you regenerate the features you will need to finish the process, since the regeneration will overwrite the currently present, already ready to use feature lists.

The final files containing all the features can be found under `ligs/csd_pool.csv` for the CSD ligands and under `ligs/lit_pool.csv` for the literature ligands. 

### Run whole model selection pipeline:
```bash
# possible modes: 0 -> oa, 1- > cp, 2 -> cc, 3 -> da_f
python main_full_model.py 0
```
This process takes around 3.5 minutes on an Intel® Core™ i7-9700K CPU. 

### Screen for candidate ligands for the OA reaction:
```bash
# possible modes: 1 -> csd ligands, 2 -> literature ligands  
python main_pool_cand.py 1
```
Besides scores of the chosen model, a list of ligands sorted by decreasing Expected Improvement (EI) values is obtained.   
