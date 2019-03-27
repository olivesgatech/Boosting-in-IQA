# Boosting-in-IQA
Code for paper - D. Temel and G. AlRegib, "Boosting in image quality assessment," 2016 IEEE 18th International Workshop on Multimedia Signal Processing (MMSP), Montreal, QC, 2016, pp. 1-6.


<p align="center">
  <img src=/Images/boosting_iqa.PNG/>
</p> 

### Paper
ArXiv: https://arxiv.org/abs/1811.08429

IEEE: https://ieeexplore.ieee.org/document/7813335

This is a brief explanation and demonstration of boosting in image quality assessment.


### Citation
If you find our paper and repository useful, please consider citing our paper:  
```
@INPROCEEDINGS{7813335, 
author={D. {Temel} and G. {AlRegib}}, 
booktitle={2016 IEEE 18th International Workshop on Multimedia Signal Processing (MMSP)}, 
title={Boosting in image quality assessment}, 
year={2016}, 
pages={1-6}, 
doi={10.1109/MMSP.2016.7813335}, 
ISSN={2473-3628}, 
month={Sep.},}

```
### Code
#### Run:
* Table 3: Run 'mslBoostingNN' 
* Table 4: Run 'mslBoostingSVR'
* Table 2: After you run the scripts above and make sure the mat files mentioned below are in your workspace ->run 'combine_NN_SVR'
* Mat files: "corr_L, corr_M, corr_T, corr_L_SVR, corr_M_SVR, corr_T_SVR, mse_L, mse_M, mse_T,
mse_L_SVR, mse_M_SVR, mse_T_SVR, P_corr_L, P_corr_M, P_corr_T, P_corr_L_SVR, P_corr_M_SVR, P_corr_T_SVR" 
#### Folder structure:
* data_L: contains the quality estiamtes of existing image quality assessment algorithms for the LIVE database
* data_M: contains the quality estiamtes of existing image quality assessment algorithms for the MULTI database
* data_T: contains the quality estiamtes of existing image quality assessment algorithms for the TID database
* mat: contains the subjective scores in the LIVE (dmos_live), the MUTLI (Multi_GT), and the TID (TID_GT) databases, the quality estiamtes of the UNIQUE method (unique_live, unique_multi, unique_tid), all the data that is used while reporting the results in the paper (NN11_Part1_3DB, SVR_Part1_3DB, NN11_Part2_3DB, SVR_Part2_3DB), orders of the methods in terms of performance (orders), table_nn (data from table 3), table_svr (data from table 4)
#### Code files:
* mslBoostingNN: Neural Netowrk experiments 
* mslBoostingSVR: Suport Vector Regression Experiments
* mse_1D: measures Mean Square Erro between two arrays
* regressMethods: regresses quality estiamtes with respect to ground truth
8 combine_NN_SVR: combines the results of NN and SVR to obtain the results in Table 2
#### Required MATLAB packages
* Neural network
* Bioinformatics ->Microarray analysis ->Expression Analysis
* Statistics and Machine Learning -> Regression -> Support Vector Machine Regression


### Abstract 
In this paper, we analyze the effect of boosting in image quality assessment through multi-method fusion. On the contrary of existing studies that propose a single quality estimator, we investigate the generalizability of multi-method fusion as a framework. In addition to support vector machines that are commonly used in the multi-method fusion studies, we propose using neural networks in the boosting. To span different types of image quality assessment algorithms, we use quality estimators based on fidelity, perceptually-extended fidelity, structural similarity, spectral similarity, color, and learning. In the experiments, we perform k-fold cross validation using the LIVE, the multiply distorted LIVE, and the TID 2013 databases and the performance of image quality assessment algorithms are measured via accuracy-, linearity-, and ranking-based metrics. Based on the experiments, we show that boosting methods generally improve the performance of image quality assessment and the level of improvement depends on the type of the boosting algorithm. Our experimental results also indicate that boosting the worst performing quality estimator with two or more methods lead to statistically significant performance enhancements independent of the boosting technique and neural network-based boosting outperforms support vector machine-based boosting when two or more methods are fused.



### Contact:

Ghassan AlRegib:  alregib@gatech.edu, https://ghassanalregib.com/, 

Dogancan Temel: dcantemel@gmail.com, http://cantemel.com/


