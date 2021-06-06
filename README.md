# Dstortion-free PCA
Distortion-free PCA on Sample Space for Highly Variable Gene Detection from Single-cell RNA-seq Data

## Usage
* demo.m : main script
* Find_HVG.m : function of the proposed method
* run_umap.m : function of the UMAP

### Parameters
* param.dim : The number of dimensions to use selecting genes (default = 3)
* para.s    : The number of selecting genes (default = 200)
* test      : Normality test for principal components
  * 'anderson' : Chi-Square goodness-of-fit test
    * param.alpha : Significance level of hypothesis test (default = 1e-4)
  * 'ske-kur'  : Skewness and kurtosis test
    * param.skewness    : Selecting the principal component in the absolute range of that parameter (default = 0.5)
    * para.kurtosis_min : Selecting the principal component that the kurtosis value larger than that parameter (default = 3)
  
### Test dataset
Test dataset is set to the sample-feature matrix (simulation_dataset.mat) with the label information (sample_labels.mat).
