library(e1071)
library(plyr)
library(reshape)
##-----------------------------------------------------------------------------------------------------------
#	           <Main functions to compute VSCS/LBMs and VSCS-AS algorithm>
##-----------------------------------------------------------------------------------------------------------

##-------------------------------------------------------------
#  Funcyion 1. This function computes the VSCS(for 90,95,99 conf. level)
#  Note: the examples used here are GLMs but can be easily adapted to
#        other models    
##-------------------------------------------------------------
ecs <- function(full.data, full.family)
{
  #full.dta: dataframe of y and x
  #full.family: family of the data, e.g. "binomial", "gaussian", "poisson"...

  #generate exhaustive set of all predictors, 2^p
  all.model<-bincombinations(p)[-1,]

  full.model <- names(full.data[,-1])
  #the full model
  full.glm <- glm(y ~0+ ., family=full.family, data = full.data)
  test.stat<-rep(NA,2^p-1)
  beta.mat<-matrix(-1,nrow=2^p-1,ncol=p)
  for(i in 1:(2^p-1)) {
   model.i<-full.model[which(all.model[i,]!=0)]
   single.glm<-glm(y ~0+ ., family = full.family, data = full.data[c(model.i, "y")])
   test.stat[i]<-anova(single.glm, full.glm, test = "Chi")$"Pr(>Chi)"[2]
  }

  ind <- which(is.na(test.stat))
  #assign 1 to these models
  test.stat[ind] <- 1

  survived.models.60 <- which(test.stat > 0.4)
  survived.models.65 <- which(test.stat > 0.35)
  survived.models.70 <- which(test.stat > 0.3)
  survived.models.75 <- which(test.stat > 0.25)
  survived.models.80 <- which(test.stat > 0.2)
  survived.models.85 <- which(test.stat > 0.15)
  survived.models.90 <- which(test.stat > 0.1)
  survived.models.95 <- which(test.stat > 0.05)
  survived.models.99 <- which(test.stat > 0.01)

  # set of survived models
  VSCS.60 <- all.model[survived.models.60,]
  VSCS.65 <- all.model[survived.models.65,]
  VSCS.70 <- all.model[survived.models.70,]
  VSCS.75 <- all.model[survived.models.75,]
  VSCS.80 <- all.model[survived.models.80,]
  VSCS.85 <- all.model[survived.models.85,]
  VSCS.90 <- all.model[survived.models.90,]
  VSCS.95 <- all.model[survived.models.95,]
  VSCS.99 <- all.model[survived.models.99,]
  return(list(p.val=test.stat, VSCS.60=VSCS.60,VSCS.65=VSCS.65,VSCS.70=VSCS.70,VSCS.75=VSCS.75,VSCS.80=VSCS.80,VSCS.85=VSCS.85,VSCS.90=VSCS.90,VSCS.95=VSCS.95,VSCS.99=VSCS.99))
}



#--------------------------------------------------------------------
## Function 2: This function computes the LBMs given the VSCS 
#--------------------------------------------------------------------

lbm <- function(VSCS)
{
p=dim(VSCS)[2]
LBM=NULL
for(i in dim(VSCS)[1]:1) {
matches <- apply(VSCS, 1, function(x) (sum(VSCS[i,]-x<0))==0)
         if(length(which(matches == TRUE)) == 1) {
           LBM<-matrix(rbind(LBM, VSCS[i,]),ncol=p)
                                                 }
                          }
return(LBM)
}



#-----------------------------------------------------------------
#   Function 3: This function generate the simulated GLMs data
#-----------------------------------------------------------------

glm.data <- function(n, p, coef.0, rho, family){
#this function returns the data as a data frame
#specify number of n samples, p predictors, rho correlation, 
#model_setup, family (as string i.e family = "binomial")
#beta0: p*1 vector of coefficients
                                                     
#covariance matrix:  Toeplitz 
  sigma <- matrix(rep(0, p*p), p, p)
  #generates the covariance matrix according to inputs
    for(i in 1:p){
      for(j in 1:p){
        sigma[i,j] <- rho^(abs(i-j))
      }
    }
   
  
  #Generate x from gaussian dist. change here if you want other form!
  x <- mvrnorm(n, rep(0, p), sigma)
  
  #initial stock response(gaussian)
  z <- x %*% as.matrix(coef.0)
  
  x=as.data.frame(x)
  
  glm.y <- function(fam, z, n, x){
  #specify the family
  #compute the response w.r.t family link
  if(fam == "binomial") #binomial dist
    {
       pr <- 1/(1+exp(-z))
       y <- rbinom(n,size=1,prob=pr) 
       return(as.data.frame(cbind(y, x)))
    }
  else if(fam == "gaussian") #gaussian dist
    {
      cov.mat <- diag(1, n, n)
      y <- mvrnorm(1, z, cov.mat)
      return(as.data.frame(cbind(y, x)))
    #normal placeholder
    }
  else if(fam == "poisson") #poisson dist
  {
     y <- rpois(n,lambda=exp(z))
     return(as.data.frame(cbind(y, x)))
  }  
  else if(fam == "Gamma") #Gamma dist
  {
     y <- rgamma(n,shape=rep(1,n),scale=1/exp(z))
     return(as.data.frame(cbind(y, x)))
  }
  else {
  return("You have entered an invalid family.")
  }
}

  #returns design matrix + response(as y)
  full.data <- glm.y(fam=family, z, n, x)
  
  
  #returns both the data and the true model
  return(list(full.data=full.data, true.model=as.numeric(c(coef.0!=0))))
  
}





#------------------------------------------------------------------
#   Function 3: This function implements the VSCS-AS
#------------------------------------------------------------------
im.sampling <- function(p, rho, B, len, omega, alpha.ast, data.glm){ 
#p: dimension
#rho: 
#omega: smoothing parameter
#alpha.ast: initial small alpha
#B:  gernerate B models given the weight
#len: length of the chain, could be replaced by some other stopping rule!!


#initialization
w <- matrix(0,len,p) #weight vectors, len*p matrix 
alpha.seq <- rep(0,len)
w[1,] <- w0 
p.val <- matrix(0,B,len)

ct<-1 # counter
gen.models<-NULL

while(ct < len) {
print(ct)
S <- matrix(0,B,p)
Gamma<-NULL
survival.no<-NULL
S <- sapply(w[ct,],rbinom,n=B,size=1)
p.val[,ct] <- apply(S, 1, glm.test,full.glm=data.glm) 

gen.models <- rbind(gen.models,S) 

alpha.temp <- sort(p.val[,ct])[B*(1-rho)]
alpha.seq[ct] <- min(alpha.temp, alpha.ast) 

survival.no <- which(p.val[,ct] > alpha.seq[ct])   ## alpha.seq[ct]?
Gamma <- S[survival.no,]
d <- dim(Gamma)[1] 
w[ct+1,] <- apply(Gamma,2,sum)/d
w[ct+1,] <- omega*w[ct+1,]+(1-omega)*w[ct,]
ct <- ct+1 
#proc.time()
}

S <- sapply(w[len,],rbinom,n=B,size=1)
p.val[,len] <- apply(S, 1, glm.test,full.glm=data.glm)
gen.models <- rbind(gen.models,S)


return(list(w=w, alpha.seq=alpha.seq, p.val=p.val,gen.models=gen.models))
}



#-----------------------------------------------------------------
# function 5: this function compute the p-val for a given model
# Note: this function will be used in function im.sampling
#-----------------------------------------------------------------
glm.test <- function(model.num, full.glm){
p1 = length(model.num)
full.variable <- names(full.glm$coefficients)[-1]  
#Note: delete "[-1]“ if GLMs does not have a intercept term 
current.model <- full.variable[which(model.num!=0)]
if (sum(model.num)==p1 ){
p.val=1
}else{
glm.current =glm(y ~ ., family =full.glm$family, data = full.data[c("y", current.model)])
p.val=anova(glm.current, full.glm$family, test= "Chi")$"Pr(>Chi)"[2]
}
return(p.val)
}

