## Connectome Based Predictive Modeling scripted in R

Code based on previously validated Python code created by Dr. Fengdan Ye. Scripts were redesigned for greater flexbility in model building, prediction, and extraction of networks. Code was NOT rigorously tested across different datasets to identify potential errors during the run.  

Currently, all scripts have parameter settings that are set to a Shen Atlas Parcellation scheme(https://pubmed.ncbi.nlm.nih.gov/23747961/) but can be updated for any atlas. This would require files similar to those contained in the ShenAtlasConfig directory here. Repeated cross-vaidation with permutation testing is available for logistic and linear models and includes  the ability to add covariates with partial correlations. 


