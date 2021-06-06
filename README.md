# Dstortion-free PCA
Distortion-free PCA on Sample Space for Highly Variable Gene Detection from Single-cell RNA-seq Data![image](https://user-images.githubusercontent.com/17615872/120912195-38954e80-c6c8-11eb-9852-616d687c5843.png)

## Usage
### Parameters
* param.dim : The number of dimensions to use selecting genes (default = 3)
* para.s    : The number of selecting genes (default = 200)
* test      : Normality test for principal components
  * 'anderson' : Chi-Square goodness-of-fit test
    * param.alpha : Significance level of hypothesis test (default = 1e-4)
  * 'ske-kur'  : Skewness and kurtosis test
    * param.skewness    : Selecting the principal component in the absolute range of that parameter (default = 0.5)
    * para.kurtosis_min : Selecting the principal component that the kurtosis value larger than that parameter (default = 3)
  
