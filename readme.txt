README FILE FOR THE PAPER 'MODEL CONFIDENCE BOUNDS FOR VARIABLE SELECTION'

1. All the figures in the paper and supplementary can be generated using these R code. And at the start of each R code file, it has the instruction about which figure does it generates.

2. At the start of each R code file, there is 'parameter setting' part, which is used to set the simulation setting. And the meaning for each parameter in the parameter setting is given below:

size: sample size 
p: the number of total covariates in the model
true.mean: the mean value of the noise 
decay.factor: the decay factor of the coefficient in the linear model
sd.epsi: the standard deviation of the noise
rho: the correlation between corvariate
r: the number of bootstrap samples
num_true_var: the number of true parameter
var_dis: the distribution of the random noise
model_type: the type of corvariates structure setting (model_type=1: means power decay of the corvariates coefficient, model_type=2: means corvariates are divided into 2 blocks, each block size is p/2 and between each block, there is no correlation between corvariate)

3. To running each R code file, first decide the simulation setting and then run the code.
