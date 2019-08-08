#### Whiten and Centering Functions ####
# Based on function constructed on matlab by  Zhilin Zhang March 13 2006 
# "Morpholically Constrained ICA for extracting weak temporally correlated signals" Neurocomputing

#1. Nargin function to check number of inputs passed to a function
nargin<-function(...){ 
  
  print(match.call()) 
  nargin <- length(as.list(match.call())) -1 
  print(nargin) 
} 



#2. Covariance Matrix Function
cov_matrix<-function(X){
  m=ncol(X) #Number of columns of input matrix
  n=nrow(X) #Number of rows of input matrix
  
  #Compute X-E(X) for each element of the matrix
  m_mean=matrix(data=1, nrow=nrow(X))%*%colMeans(X)
  D=X-m_mean
  
  #Computing each element of the matrix
  C<-(n)^-1*t(D)%*%D #Equivalent to computing [X-E(X)]^2
  return(C)
}


#3. Centering
centering_ica<-function(X){
  m_mean=matrix(data=1, nrow=nrow(X))%*%colMeans(X)
  D=X-m_mean
  return(D)
}


#4. Whiten matrix based on Hyvarinen & Oja whitening in fastICA
whiten_ica2<-function(X){
  #Centering
  X_c=centering_ica(X)
  
  n=nrow(X_c)
  t=ncol(X_c)
  #Estimate E and D
  epsilon = 1e-13 #To prevent NaNs
  E=eigen(cov_matrix(X_c))$vectors
  D=eigen(cov_matrix(X_c))$values+epsilon
  D_sqrt=D^(-0.5)

  #Whitening matrix
  W=E%*%diag(D_sqrt)%*%t(E)
  #Transformed input
  X_hat=X_c%*%W
  
  output=list("E"=E, "D"=D ,"D_sqrt"=D_sqrt ,"W"=W, "X_hat"=X_hat)
  return (output)
  }

#5. Lu Rajapakse cICA Estimation
Lu_Rajapakse_cICA<-function(X, ref, threshold, w0, learningRate, mu0, lambda0, gamma, maxIter, OverValue, verbose=TRUE){

#Define dimensions of ICs
ICnum=dim(X)[2]
IClen=dim(X)[1]

#Define starting weight vectors (y=w'X inverse matrix)
w=w0
oldw=w
#Define Lagrange Multipliers for inequality constraint (g(y:W))
mu=mu0
#Define Lagrange Multipliers (h(y:W))
lambda=lambda0
#Define flag and loop parameters
flag=1
loop=1

#Compute Covariance Matrix
Rxx=cov_matrix(X)

#While loop
while (flag==1) {
  
  #Output using current iteration
  y=X%*%w
  
  ###Calculate the first order deviation of the Lagarange function
  std_y=sd(y) #Compute Standard Deviation of y
  v_gauss=rnorm(IClen, mean=0, sd=std_y) #Compute Gaussian rv with same mean and variance as X
  #rou parameter of the optimization problem
  rou=mean(log(cosh(y))-log(cosh(v_gauss)))
  #First order derivative of the augmented Lagrange Function
  LL=sign(rou) * ((t(tanh(y))%*%X)/IClen) - mu*(t(y-ref)%*%X)/IClen - lambda*(t(y)%*%X)/IClen
  
  ###Calculate second order deviation of Lagrange function
  Sw = sign(rou) * mean(1-tanh(y)^2) - mu - lambda
  
  #update weight vector
  w=w - (t(learningRate*(LL%*%solve(Rxx))))/Sw
  #Normalize weight vector
  w=w/norm(X, type=c("F"))
  
  #Update Mu
  thr=threshold*(1-exp(-loop))
  #Standarize y and ref signal as required 
  y_stand=(y-mean(y))/sd(y)
  ref_stand=(ref-mean(ref))/sd(ref)
  #Inequality Constraint //TRY DIFFERENT DISTANCE FUNCTIONS
  g=mean((y_stand-ref_stand)^2)-thr
  #g=diss.DTWARP(y_stand, ref_stand)-thr #Time Warping Distance
  #g=1/(abs(mean(y_stand-ref_stand)))^2 - thr
  mu=max(0, mu + gamma * g)
  #mu=max(gamma*g, -mu/gamma)
  #Equality Constraint
  h=mean((y^2)-1)
  lambda = lambda + gamma*h
  
  #Convergence Criterion
  wchange = 1-as.numeric(abs(t(oldw)%*%w))
  if(verbose){print(sprintf("No.%d iteration: change in w is %g",loop, wchange))}
  
  if(wchange<OverValue){
    if(verbose){print(sprintf("Converged after %d iteration",loop))}
    flag=0
  }
  
  if(loop>=maxIter){
    if(verbose){print(sprintf("After %d iteration, still cannot converge",loop))}
    flag=0
  }
  
  oldw = w
  loop = loop + 1
  
} 
#End of While loop
y=X%*%w
if(verbose){print("End of cICA algorithm")}

output=list(y,w)
return(output) 

}

#6. SNR (Signal-to-Noise ratio)
SNR<-function(y,s){
  
  MSE=mean((y-s)^2)
  vsource=var(s)^2
  snr=10*log10(vsource/MSE)
  
  return (snr)
  }


#Test matrixes
#S <- matrix(runif(1000), 500, 2) 
#A <- matrix(c(0.31, 0.12, -1.34, 3.25), 2, 2, byrow = TRUE) 
#X <- S %*% A
#library("entropy", lib.loc="~/R/win-library/3.5")

#Test matrixes Simulation 1: Random ref signals
N=300
k=seq(1:N)

ts = 1e-4 # sampling period
f1 = 0.061/ts    # true frequency
f2 = 0.054/ts    
f3 = 0.028/ts

S=matrix(0, N, 5)


S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k))
S[,2] = cos(2*pi*f2*ts*k)
S[,3] = cos(2*pi*f3*ts*k + 2)
S[,4] = runif(N)
S[,5] = runif(N)

S=whiten_ica2(S)$X_hat
A = matrix(runif(25), 5,5)
X = S%*%A

#Set Parameters for cICA
w = as.matrix(runif(5, min=0, max=1))
w=w/norm(w, type=c("F"))

mu0 = 1
lambda0 = 1
gamma = 1
learningRate = 1
OverValue=0.000001  
maxIter = 200000
threshold = 1.75

ref1 = 2*sin(2*pi*f1*ts*k)
ref2 = runif(N)

c1=cICA(whiten_ica2(X)$X_hat, ref1, threshold, w, learningRate, mu0, lambda0, gamma, maxIter, OverValue)

#Test matrixes Simulation 2: Clustered ref signals
library("TSclust", lib.loc="~/R/win-library/3.5")
library("factoextra", lib.loc="~/R/win-library/3.5")
library("tidyverse", lib.loc="~/R/win-library/3.5")
library("dendextend", lib.loc="~/R/win-library/3.5")
#1. Apply whitening to data matrix
X_whiten<-whiten_ica2(X)$X_hat

#2. Construct distance matrix using Integrated Periodgram Distance
Dist<-diss(t(X), "INT.PER")
fviz_dist(Dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#3. Evaluate different Hierarchical clustering algorithms (Agglomerative)
m<-c( "average", "single", "complete", "ward")
names(m)<-c( "average", "single", "complete", "ward")
ac <- function(x) {agnes(Dist, method = x)$ac}
map_dbl(m, ac) #Ward's yields the best results

#4. Visualize Dendrogram and generate subgroups
hc<-agnes(Dist, method = "complete")
pltree(hc, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
clusters<-cutree(as.hclust(hc), k = 3)
table(clusters)
#Visualize cluster
fviz_cluster(list(data = Dist, cluster = clusters))

#5. Find Centroids
centroids<-as.data.frame(t(apply(t(X), 2, function (x) tapply(x, clusters, mean))))

#6. Define New Reference as Centroid
ref3=as.numeric(centroids[,1])

#7. Run cICA
c2=cICA(whiten_ica2(X)$X_hat, ref3, threshold, w, learningRate, mu0, lambda0, gamma, maxIter, OverValue)

#Check against fastICA
library("fastICA", lib.loc="~/R/win-library/3.5")
ic=fastICA(X, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)

#Check against PCA
library("FactoMineR", lib.loc="~/R/win-library/3.5")
pc=prcomp(X, scale=TRUE)

#Compute Normalize RMSE
library("hydroGOF", lib.loc="~/R/win-library/3.5")
nrmse(as.numeric(unlist(c2[1])), as.numeric(S[,1]), norm="maxmin")
nrmse(as.numeric(unlist(c1[1])), as.numeric(S[,1]), norm="maxmin")
nrmse(as.numeric(ic$S[,1]), as.numeric(S[,1]), norm="maxmin")
nrmse(as.numeric(pc$x[,1]), as.numeric(S[,1]), norm="maxmin")
