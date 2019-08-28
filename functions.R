##################################
# MAIN FUNCTIONS FOR THE 'CONFIDENCE BOUNDS FOR VARIABLE SELECTION'
##################################
library(parallel)
library(methods)
library(leaps)
library(lars)
library(MASS)
library(glmnet)
library(ncvreg)
library(parcor)
library(lqa)


#Generate data function
Generate.Data<-function(size,p, true.mean, true.sd, decay.factor,sd.epsi,rho,num_true_var, model_type, var_dis)
{
  #multinominal correlated simulation
  omiga<-matrix(0,nrow=p,ncol=p)
  if (model_type == 1){
      for(i in 1:p)
      {
        for(j in 1:p)
        {
          omiga[i,j]<-rho^(abs(i-j))
        }
      }
  }else if(model_type == 2){ # setting2 generate data with 2 block, each block size is p/2 and within each block, there is rho correlation
      for(i in 1:p){
        for(j in 1:p)
        {
          if (i == j) omiga[i,j] <- 1 else{
            if (i <= p/2 && j <= p/2) {omiga[i,j]<-rho}
            if (i > p/2 && j > p/2) {omiga[i,j]<-rho}
          }

        }
      }
  }

  data<-mvrnorm(n=size,mu=rep(0,p),omiga)
  colnames(data)=paste("x",1:p,sep="")
  ############################ ERROR USING T distribution
  if (var_dis == 't') u<-rt(size, df)
  ###########################ERROR USING NORMAL DISTRIBUTION
  if (var_dis == 'normal') u<-rnorm(size,mean=0,sd = sd.epsi)
  ###########################ERROR USING CAUCHY DISTRIBUTION
  if (var_dis == 'cauchy') u<-rcauchy(size,location=0,scale=1)
  # Laplase distribution
  if (var_dis == 'Laplase') u<-rdoublex(size,mu=0,lambda =sqrt(sd.epsi^2/2))
  b <- vector(length = num_true_var, mode = "numeric")
  for(i in 1:num_true_var){
    b[i] = 1 * decay.factor^(i-1)
  }
  y<- data[,1:num_true_var] %*% b + u
  return(data.frame(data,y))
}

#Transfrom a list of variable name denoted selection results to 0-1 matrix result  
f01<-function(object){  
  full<-paste("x",1:p,sep="")
  matrix.01<-data.frame(matrix(0,length(object),p));colnames(matrix.01)<-full
  for(i in 1:length(object)){matrix.01[i,(full %in% object[[i]])]<-1}
  return(matrix.01)
}

#Transform 0-1 matrix variable selection result ot a list of variable name denoted selection resuls
f02<-function(object) 
{
  full<-paste("x",1:p,sep="")
  var.matrix<-vector(mode="list",length=dim(object)[1])
  for(i in 1:dim(object)[1])
  {
    var.matrix[[i]]<-full[object[i,]>0]
  }
  return(var.matrix)
}

#big_matrix_illustration function:numerical vector to name list transformation
f03 <- function(object){
  full <- paste("x",1:p,sep="")
  name_vector <- full[object > 0]
}

#function3
inter.union<-function(object){
  j=3;
  inter.sect<-intersect(object[[1]],object[[2]])
  union.sect<-union(object[[1]],object[[2]])
  while(j<=k){
    inter.sect<-intersect(inter.sect,object[[j]])
    union.sect<-union(union.sect,object[[j]])
    j<-j+1
  } 
  return(list(inter<-inter.sect,union<-union.sect))  
}

#online adaptive lasso code:http://www4.stat.ncsu.edu/~boos/var.select/lasso.adaptive.html
lasso.adapt.bic2<-function(x,y){
  
  # adaptive lasso from lars with BIC stopping rule 
  # this one uses the "known variance" version of BIC with RSS/(full model mse)
  # must use a recent version of R so that normalize=FALSE can be used in lars
  
  require(lars)
  ok<-complete.cases(x,y)
  x<-x[ok,]                            # get rid of na's
  y<-y[ok]                             # since regsubsets can't handle na's
  m<-ncol(x)
  n<-nrow(x)
  x<-as.matrix(x)                      # in case x is not a matrix
  
  #  standardize variables like lars does 
  one <- rep(1, n)
  meanx <- drop(one %*% x)/n
  xc <- scale(x, meanx, FALSE)         # first subtracts mean
  normx <- sqrt(drop(one %*% (xc^2)))
  names(normx) <- NULL
  xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
  
  out.ls=lm(y~xs)                      # ols fit on standardized
  beta.ols=out.ls$coeff[2:(m+1)]       # ols except for intercept
  w=abs(beta.ols)                      # weights for adaptive lasso
  xs=scale(xs,center=FALSE,scale=1/w)  # xs times the weights
  object=lars(xs,y,type="lasso",normalize=FALSE)
  
  # get min BIC
  # bic=log(n)*object$df+n*log(as.vector(object$RSS)/n)   # rss/n version
  sig2f=summary(out.ls)$sigma^2        # full model mse
  bic2=log(n)*object$df+as.vector(object$RSS)/sig2f       # Cp version
  step.bic2=which.min(bic2)            # step with min BIC
  
  fit=predict.lars(object,xs,s=step.bic2,type="fit",mode="step")$fit
  coeff=predict.lars(object,xs,s=step.bic2,type="coef",mode="step")$coefficients
  coeff=coeff*w/normx                  # get back in right scale
  st=sum(coeff !=0)                    # number nonzero
  mse=sum((y-fit)^2)/(n-st-1)          # 1 for the intercept
  
  # this next line just finds the variable id of coeff. not equal 0
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind<-0
  intercept=as.numeric(mean(y)-meanx%*%coeff)
  return(list(fit=fit,st=st,mse=mse,x.ind=x.ind,coeff=coeff,intercept=intercept,object=object,
              bic2=bic2,step.bic2=step.bic2))
}

#Get the bootsrap variable selection result
# Note you can use this function for different variable selection method, but the default one is adaptive lasso
BOOT.CI<-function(x, dep.var.index, r){
  var.instances<-vector(mode="list",length=r)
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    boot.data<-x[ind,]
    # fit<-regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
    # var.instances[[j]]<-names(coef(fit,which.min(summary(fit)$bic)))[-1]
    fit<-lasso.adapt.bic2(x=boot.data[,1:p],y=boot.data$y)
    var.instances[[j]]<-names(fit$coeff)[(fit$coeff)!=0]
    # adalasso<-adalasso(X=as.matrix(boot.data[,1:p]),y=boot.data$y,k=10)
    # var.instances[[j]]<-full.var[adalasso$coefficients.adalasso!=0]
    # opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    # lasso.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    # beta4<-coef(lasso.fit,s=(opt.lambda+0.2))[,1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta4)[beta4!=0]]
    # fit <- cv.lqa(y.train = boot.data[,dep.var.index],x.train = as.matrix(boot.data[,-dep.var.index]),intercept = TRUE,lambda.candidates = list(c(1)),family = gaussian(),penalty.family = lasso,n.fold = 10,loss.func = "aic.loss")
    # beta<- coef(fit$best.obj)[-1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var.instances)
}

# GET the residual bootsrap variable selection results by using adaptive lasso and 
# you can give a specified lambda
RES.BOOT.CI2<-function(x, dep.var.index, r,lmbd){
  var.instances<-vector(mode="list",length=r)
  ## = adaLASSO = ##
  tau=1
  lasso_init=glmnet(as.matrix(x[,1:p]),x$y)
  first.step.coef=lasso_init$beta[,which.min(abs(lasso_init$lambda-lmbd))]
  penalty.factor=abs(first.step.coef+1/sqrt(nrow(x)))^(-tau)
  adalasso=glmnet(as.matrix(x[,1:p]),x$y,penalty.factor=penalty.factor)
  beta_est=adalasso$beta[,which.min(abs(adalasso$lambda-lmbd))]
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta_est
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta_est
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    lasso_init_boot=glmnet(as.matrix(boot.data[,1:p]),boot.data$y)
    first.step.coef_boot=lasso_init_boot$beta[,which.min(abs(lasso_init_boot$lambda-lmbd))]
    penalty.factor_boot=abs(first.step.coef_boot+1/sqrt(nrow(x)))^(-tau)
    adalasso_boot=glmnet(as.matrix(boot.data[,1:p]),boot.data$y,penalty.factor=penalty.factor_boot)
    beta_est_boot=adalasso_boot$beta[,which.min(abs(adalasso_boot$lambda-lmbd))]
    var.instances[[j]]<-full.var[beta_est_boot!=0]
  }
  return(var.instances)
}

#Get the residual bootsrap variable selection result by using adaptive lasso
# Note you can use this function for different variable selection method, but the default one is adaptive lasso
RES.BOOT.CI<-function(x, dep.var.index, r){
  var.instances<-vector(mode="list",length=r)
  adalasso<-adalasso(X=as.matrix(x[,1:p]),y=x$y,k=10)
  res_original <- x$y - as.matrix(x[,1:p]) %*% adalasso$coefficients.adalasso
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% adalasso$coefficients.adalasso
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    # fit<-regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
    # var.instances[[j]]<-names(coef(fit,which.min(summary(fit)$bic)))[-1]
    # fit<-lasso.adapt.bic2(x=boot.data[,1:p],y=boot.data$y)
    # var.instances[[j]]<-names(fit$coeff)[(fit$coeff)!=0]
    adalasso_tem<-adalasso(X=as.matrix(boot.data[,1:p]),y=boot.data$y,k=10)
    var.instances[[j]]<-full.var[adalasso_tem$coefficients.adalasso!=0]
    # opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    # lasso.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    # beta4<-coef(lasso.fit,s=(opt.lambda+0.2))[,1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta4)[beta4!=0]]
    # fit <- cv.lqa(y.train = boot.data[,dep.var.index],x.train = as.matrix(boot.data[,-dep.var.index]),intercept = TRUE,lambda.candidates = list(c(1)),family = gaussian(),penalty.family = lasso,n.fold = 10,loss.func = "aic.loss")
    # beta<- coef(fit$best.obj)[-1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var.instances)
}


# get residual bootstrap variable selection models from SCAD and MCP
RES.BOOT.CI3 <- function(x, dep.var.index, r, pnlt){
  var<-vector(mode="list",length=r)
  fit<-cv.ncvreg(X=x[,1:p],y=x$y,penalty=pnlt)
  fit2<-fit$fit
  beta<-fit2$beta[,fit$min]
  beta<-beta[2:(p+1)]
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data$y,penalty=pnlt)
    fit2_tem<-tem$fit
    beta_tem<-fit2_tem$beta[,tem$min]
    beta_tem <- beta_tem[2:(p+1)]
    var[[j]]<-names(beta_tem)[which(beta_tem!=0)]
  }
  return(var)
}

# Use LAD and Sqrt method to get residual bootstrap variable models
RES.BOOT.CI4 <- function(x, dep.var.index, r, q, beta_index){
  var<-vector(mode="list",length=r)
  fit<-slim(X=as.matrix(x[,1:p]),Y=x$y,method="lq",q=q)
  beta <- fit$beta[,beta_index]
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- slim(X=as.matrix(boot.data[,1:p]),Y=boot.data$y,method = 'lq', q = q)
    var[[j]]<-full.var[tem$beta[,beta_index]!=0]
  }
  return(var)
}

# GET residual boostrap variable selection models with stepwise BIC
RES.BOOT.CI5 <- function(x, p, r){
  var<-vector(mode="list",length=r)
  fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
  if (names(coef(fit,which.min(summary(fit)$bic)))[1] == '(Intercept)'){
    beta <- vector(mode = 'numeric', length = p)
    beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))[-1]
  } else {
    beta <- vector(mode = 'numeric', length = p)
    beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))
  }
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  for(j in 1:r) { 
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
    if (names(coef(tem,which.min(summary(tem)$bic)))[1] == '(Intercept)'){
      var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))[-1]
    } else {
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))
    }
    
  }
  return(var)
}

# Get boostrap models from the modified lasso
BOOT.MODI.LASSO <- function(x, dep.var.index, r, threshold){
  var_lasso<-vector(mode="list",length=r)
  for (i in 1:r){
    ind=sample(1:nrow(x),nrow(x),replace=T)
    boot.data<-x[ind,]
    opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    lasso.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    beta<-coef(lasso.fit,s=opt.lambda)[,1][which(abs(coef(lasso.fit,s=opt.lambda)[,1]) > threshold)]
    var_lasso[[i]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var_lasso)
}


#Calculate the Model confidences bounds and their corresponding coverage rate
CI<-function(var.list,var.matrix,p)
{
  full.var<-colnames(var.matrix)
  colsum<-apply(var.matrix,2,sum)
  order<-order(colsum,decreasing = T)
  freq<-vector(length=p+1);freq[1]<-0
  lower<-vector(mode="list",length=p+1)
  upper<-vector(mode="list",length=p+1)
  for(i in 0:p)
  { 
    cap<-vector(length=p+1);cap[1]<-0  
    for(j in 0:(p-i))
    {
      if (j==0 & i!=0)
      {
        uppertest<-full.var[order[1:i]]
        for(r in 1:length(var.list))
        {
          if(all(var.list[[r]]%in%uppertest)) cap[j+1]<-cap[j+1]+1
        }
      }else{
        if(j!=0){
          lowtest<-full.var[order[1:j]]
          uppertest<-full.var[order[1:(j+i)]]
          for(r in 1:length(var.list))
          {
            if(all(all(lowtest%in%var.list[[r]]),all(var.list[[r]]%in%uppertest))) cap[j+1]<-cap[j+1]+1
          }
        }
        if (j == 0){
          for (r in 1:length(var.list)){
            if (identical(var.list[[r]],character(0))) cap[j+1] <- cap[j+1] + 1
          }
        }
      }
    }  
    freq[i+1]<-max(cap)/r
    maxlocation<-which.max(cap)
    if(maxlocation==1)
    {
      if (i != 0){
        lower[[i+1]]<-''
        upper[[i+1]]<-full.var[order[1:i]]
      } else if (i == 0){
        lower[[1]]<-''
        upper[[1]]<-''
      }
    }else{
      lower[[i+1]]<-full.var[order[1:(maxlocation-1)]]
      upper[[i+1]]<-full.var[order[1:(maxlocation-1+i)]]
    }
  }   
  result<-list(freq=freq,lower=lower,upper=upper)
  return(result)
}

# A function used in calculating optimal MCB
var.f<-function(k){
  full.combine<-full.var
  var.include<-vector("list",length=choose(p,k))
  for(i in 1:length(var.include)){
    var.include[[i]]<-full.combine[t(combn(p,k))[i,]]
  }
  return(var.include)
}

# A function used in calculating optimal MCB
f1<-function(varmatrix,varinstance)
{
  freq<-rep(0,length(varmatrix))
  for(i in 1:length(varmatrix))
  {
    for(j in 1:length(varinstance))
    {
      freq[i]<-freq[i]+all(varmatrix[[i]] %in% varinstance[[j]])
    }
    freq[i]<-freq[i]/length(varinstance)
  }
  return(max(freq))
}

# A function used in calculating optimal MCB
f2<-function(varmatrix,varinstance)
{
  freq<-rep(0,length(varmatrix))
  for(i in 1:length(varmatrix))
  {
    for(j in 1:length(varinstance))
    {
      freq[i]<-freq[i]+all(varinstance[[j]] %in% varmatrix[[i]])
    }
    freq[i]<-freq[i]/length(varinstance)
  }
  return(max(freq))
}

# This function is to all possiable lower bound with length x, and there is no repeat with choose((p-low.len),width)
low.var.norepeat<-function(p,low.len){ 
  var.low<-vector("list") 
  var<-as.matrix(t(combn(full.var,low.len)))
  var<-apply(var,1,paste,collapse="#")
  return(sapply(var,strsplit,split="#",USE.NAMES = F))
}

# A function used in calculating optimal MCB
low_up.var<-function(p,low.len,width){
  # This function is to all possiable lower bound with length x, and there is repeat with choose((p-low.len),width)
  low.var.repeat<-function(p,low.len,width){ 
    var.low<-vector("list") 
    var<-as.matrix(t(combn(full.var,low.len)))
    var<-apply(var,1,paste,collapse="#")
    object<-as.matrix(rep(var,each=choose((p-low.len),width)))
    return(sapply(object,strsplit,split="#",USE.NAMES = F))
  }
  low.var.object<-low.var.repeat(p,low.len,width)
  
  add.var<-function(p,low.len,width){
    var<-as.matrix(t(combn(full.var,low.len)))
    ind<-t(apply(var,1,function(var){!full.var%in%var}))
    object<-list()
    for(i in 1:nrow(ind)){
      object[[i]]<-full.var[ind[i,]]
    }  
    add<-lapply(object,combn,m=width)
    add1<-lapply(add,t);len<-nrow(add1[[1]])
    temp<-matrix(,nr=length(add1)*len,nc=ncol(add1[[1]]))  
    for(i in 1:length(add1)){
      temp[((i-1)*len+1):(i*len),]<-add1[[i]]
    }
    return(temp)
  }
  add.var.object<-add.var(p,low.len,width)
  
  up.var<-list()
  for(i in 1:length(low.var.object)){
    up.var[[i]]<-c(low.var.object[[i]],add.var.object[i,])
  }
  return(list(low=low.var.object,up=up.var))
}

# The function to calculate optimal MCB
##i is the  length of lower bound
##w is the width of the confidence bound
cioptimal<-function(var.instances, p){
  res<-data.frame(matrix(,nr=p,nc=p))
  colnames(res)<-paste("width",0:(p-1),sep="")
  for(w in 0:(p-1)){
    for(i in 1:(p-w)){
      var.set<-low_up.var(p,i,w) 
      count<-rep(0,length(var.set$low))
      cap<-vector(length=length(count))   
      for(j in 1:length(var.set$low)){
        
        for(k in 1:r){
          count[j]<-count[j]+all(all(var.set$low[[j]] %in% var.instances[[k]]),all(var.instances[[k]] %in% var.set$up[[j]]))
        }
        cap[j]<-count[j]/r     
      }
      res[i,(w+1)]<-max(cap)
    }
  }
  cap_res<-apply(res,2,max,na.rm=T)
  return(cap_res)
}

# Given the data set, to see whether MCBs could cover the true model
captrue_model<-function(x, dep.var.index, r, lambda_list)
{
  freq_matrix <- matrix(0, nrow = length(lambda_list), ncol = dep.var.index)
  truecapture_matrix <- matrix(0, nrow = length(lambda_list), ncol = dep.var.index)
  for (i in 1:length(lambda_list)){
    var.instances <- RES.BOOT.CI2(x, dep.var.index, r, lambda_list[i])
    var.01 <- f01(var.instances)
    result<-CI(var.instances,var.01,p)
    lower<-result$lower
    upper<-result$upper
    freq_matrix[i,] <- result$freq
    truecapture_matrix[i,] <- cap_true(lower,upper,p)
  }
  results<-list(freq=freq_matrix,truecapture=truecapture_matrix)
  return(results)
}

# A functional function for the function "captrue_model"
cap_true<-function(low_list,upper_list,p){
    truecapture<-vector()
    for(i in 1:(p+1))
    {
      if (low_list[[i]][1] == ""){
        if (all(all(true.model %in% upper_list[[i]]))) truecapture[i]<-1 else truecapture[i]<-0
      } else{
        if(all(all(low_list[[i]]%in%true.model),all(true.model %in% upper_list[[i]]))) truecapture[i]<-1 else truecapture[i]<-0
      } 
    }
    return(truecapture)
}

# A function a analysis_of_bootstrap_validity: used to quantify the selection method consistency rate 
judge_equal <- function(data, compared_model){
  adalasso_tem <- adalasso(X=as.matrix(data[,1:p]),y=data$y,k=10)
  var.instances <- full.var[adalasso_tem$coefficients.adalasso!=0]
  if (length(var.instances) == length(compared_model) && sum(var.instances%in%compared_model) == length(compared_model)) equal <-1 else equal <- 0
  return(equal)
}



#MODEL RANKING
model_ranking <- function(input_matrix, input_vector){
  row_dim <- dim(input_matrix)[1]
  if (row_dim == 0){
    input_matrix = rbind(input_matrix, input_vector)
  } else{
    insert_index = 1
    for (i in 1:row_dim){
      k = 0
      if (sum(input_vector) > sum(input_matrix[i,])){
        insert_index <- insert_index + 1
        k = 1
      } else if (sum(input_vector) == sum(input_matrix[i,])){
        model_dif <- input_vector - input_matrix[i,]
        for (j in 1:p){
          if (model_dif[j] == -1){
            insert_index <- insert_index + 1
            k = 1
            break
          } else if (model_dif[j] == 1) {
            k = 0
            break
          }
        }
      }
      if (k == 0){
        break
      }
    }
    if (insert_index == (row_dim + 1)){
      input_matrix <- rbind(input_matrix, input_vector)
    } else{
       input_matrix <- insertRow(input_matrix, input_vector, insert_index)
    }
  }
  return(input_matrix)
}

# A function in big_matrix_illustration.R
insertRow <- function(existingDF, newrow, r){
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}



# A revised version of MCB_compare_other, so you can set the tuning parameter of adaptive lasso by yourself
MCB_compare_other2 <- function(x, dep.var.index, r, candidate_lambda){
  if (calculate_MCB == 1){
    MCB_result <- MCB_acc_cardi_calculate(x, r, candidate_lambda)
    cardinality_MCB <- tryCatch(MCB_result$cardinality_MCB_sub, error = function(e){print('error3'); problem = MCB_result})
    capture_MCB <- tryCatch(MCB_result$capture_MCB_sub, error = function(e){print('error3'); problem = MCB_result})
    score_MCB <- tryCatch(MCB_result$score_MCB_sub,error = function(e){print('error3'); problem = MCB_result})
  }
  if (calculate_VSCS == 1){
    VSCS_result <- VSCS_acc_cardi_calculate(x)
    cardinality_VSCS <- tryCatch(VSCS_result$cardinality_VSCS_sub,error = function(e){print('error4'); problem = VSCS_result})
    capture_VSCS <- tryCatch(VSCS_result$capture_VSCS_sub,error = function(e){print('error4'); problem = VSCS_result})
  }
  if (calculate_MCB == 1 && calculate_VSCS == 1) return_result<-list(cardinality_MCB=cardinality_MCB,cardinality_VSCS=cardinality_VSCS,capture_MCB=capture_MCB,capture_VSCS=capture_VSCS,score_MCB= score_MCB)
    if (calculate_MCB == 1 && calculate_VSCS == 0) return_result<-list(cardinality_MCB=cardinality_MCB,capture_MCB=capture_MCB,score_MCB= score_MCB)
    if (calculate_MCB == 0 && calculate_VSCS == 1) return_result<-list(cardinality_VSCS=cardinality_VSCS,capture_VSCS=capture_VSCS)
    return(return_result)
}

# Calcuate the true model capture and cardinality of VSCS
VSCS_acc_cardi_calculate <- function(dataset){
  cardinality_VSCS_sub <- vector('numeric', length = length(alpha))
  capture_VSCS_sub <- vector('numeric', length = length(alpha))
  x_transform <- cbind(dataset[,p+1], dataset[,1:p])
  colnames(x_transform)[1] <- 'y'
  result_VSCS <- ecs(full.data=x_transform, full.family='gaussian')
  for (zz in 1:length(alpha)){
    capture_VSCS_sub[zz] <- VSCS_capture_true_model(result_VSCS[[zz+1]], true.model_vector)
    if (is.null(result_VSCS[[zz+1]])) cardinality_VSCS_sub[zz] <- 0 else if (is.vector(result_VSCS[[(zz+1)]])) {
      cardinality_VSCS_sub[zz] <- 1 
    } else {cardinality_VSCS_sub[zz] <- dim(result_VSCS[[zz+1]])[1]}
  }
  return_value <- list(cardinality_VSCS_sub=cardinality_VSCS_sub, capture_VSCS_sub= capture_VSCS_sub)
  return(return_value)
}

# Calculate the true model capture and cardinality and score of MCB
MCB_acc_cardi_calculate <- function(dataset, r, lambda_list){
  cardinality_MCB_sub <- matrix(0, nrow = length(lambda_list), ncol = length(alpha))
  capture_MCB_sub <- matrix(0, nrow = length(lambda_list), ncol = length(alpha))
  score_MCB_sub <- matrix(0, nrow = length(lambda_list), ncol = length(alpha))
  for (i in 1:length(lambda_list)){
    lambda <- lambda_list[i]
    var.instances<-RES.BOOT.CI2(dataset,p + 1,r,lambda)
    var.01<-f01(var.instances)
    result_MCB<-CI(var.instances,var.01,p)
    ci_lower<-result_MCB$lower
    ci_upper<-result_MCB$upper
    # print("print out MCB")
    # print(lambda_list[i])
    # print("print lower bound")
    # print(ci_lower)
    # print("upper bound")
    # print(ci_upper)
    # print(result_MCB$freq)
    for (j in 1:length(alpha)){
      talpha <- alpha[j]
      lower_alpha<-ci_lower[which(result_MCB$freq > (1 - talpha))[1]]
      upper_alpha<-ci_upper[which(result_MCB$freq > (1 - talpha))[1]]
      if (lower_alpha[[1]][1] == ""){
        if (all(all(true.model %in% upper_alpha[[1]]))) {
          capture_MCB_sub[i, j] <- 1 
          score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha)
          }else {
            capture_MCB_sub[i, j] <- 0
            score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha) + 2/talpha * (length(true.model) - length(upper_alpha))
          }
      } else{
        if(all(all(lower_alpha[[1]]%in%true.model),all(true.model %in% upper_alpha[[1]]))) {
          capture_MCB_sub[i, j] <- 1 
          score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha)
          }else {
            capture_MCB_sub[i, j] <- 0
            score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha) + 2/talpha * max((length(true.model) - length(upper_alpha)),0) + 2/talpha * max((length(lower_alpha) - length(true.model)),0)
          }
      }
      lower_alpha<-f01(lower_alpha)
      upper_alpha<-f01(upper_alpha)
      cardinality_MCB_sub[i,j] <- 2^(sum(upper_alpha) - sum(lower_alpha))
    }
  }
  return_value = list(cardinality_MCB_sub=cardinality_MCB_sub, capture_MCB_sub= capture_MCB_sub, score_MCB_sub=score_MCB_sub)
  return(return_value)
}


# A function used in MCB and VSCS comparison
VSCS_capture_true_model <- function(model_matrix, real_model_vector){
  capture <- 0
  if (is.null(dim(model_matrix))){return(capture)} else{
    if (is.vector(model_matrix)) {
      if(sum(abs(model_matrix - real_model_vector)) == 0){
        capture <- 1
        return(capture)
      } else {
        return(capture)
      }
    } else if (dim(model_matrix)[1] == 0) return(capture) else{
    for (i in 1:dim(model_matrix)[1]){
      if (sum(abs(model_matrix[i,] - real_model_vector)) == 0){
        capture <- 1
        break
        }
      }
      return(capture)
    }
  }
}

# A plot function used in comparison of different model selection method
model_selection_method_plot_function <- function(epi_dis){
  if (epi_dis == 'normal'){
    plot(c(0:p)/p,final_result[,2],lty=2,lwd=3,cex.axis=1.7, pch = 3,cex.lab=cexnum,cex.main=1.3,xlim=c(0,1),type='b',ylim=c(0,1),xlab='w/p',ylab = 'Coverage Rate')
    lines(c(0:p)/p,final_result[,1],lty=3,lwd=3)
    lines(c(0:p)/p,final_result[,3],lty=4,lwd=3)
    lines(c(0:p)/p,final_result[,4],lty=5,lwd=3)
    lines(c(0:p)/p,result_stepwise$freq,lty=1, lwd=3)
    lines(c(0:p)/p, result_LAD[[5]]$freq, type = 'b', lty =1,pch=1, lwd =3)
    lines(c(0:p)/p, result_SQRT[[5]]$freq, type = 'b',pch = 2,lty =1, lwd =3)
    legend('bottomright',legend=c("Stepwise","Lasso",'Adaptive Lasso','LAD Lasso','SQRT Lasso',"SCAD",'MCP'),lty=c(1,2,3,1,1,4,5), pch = c(NA,3,NA,1,2,NA,NA),text.width = 0.35,lwd=rep(3,7),cex = 2,pt.cex = 2)
  } else {
    plot(c(0:p)/p,final_result[,2],lty=2,lwd=3,cex.axis=1.7, pch = 3,cex.lab=cexnum,cex.main=1.3,xlim=c(0,1),type='b',ylim=c(0,1),xlab='w/p',ylab = 'Coverage Rate')
    lines(c(0:p)/p,final_result[,1],lty=3,lwd=3)
    lines(c(0:p)/p,final_result[,3],lty=4,lwd=3)
    lines(c(0:p)/p,final_result[,4],lty=5,lwd=3)
    lines(c(0:p)/p,result_stepwise$freq,lty=1, lwd=3)
    lines(c(0:p)/p, result_LAD[[2]]$freq, type = 'b', lty =1,pch=1, lwd =3)
    lines(c(0:p)/p, result_SQRT[[4]]$freq, type = 'b',pch = 2,lty =1, lwd =3)
    legend('bottomright',legend=c("Stepwise","Lasso",'Adaptive Lasso','LAD Lasso','SQRT Lasso',"SCAD",'MCP'),lty=c(1,2,3,1,1,4,5), pch = c(NA,3,NA,1,2,NA,NA),text.width = 0.35,lwd=rep(3,7),cex = 2,pt.cex = 2)
  }

}

