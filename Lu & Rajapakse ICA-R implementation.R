###Pre-Processing Functions###
library("MASS", lib.loc="C:/Program Files/R/R-3.5.1/library")
library("pracma", lib.loc="~/R/win-library/3.5")

#1. Centering
centering_icaR<-function(X){
  m_mean=matrix(data=1, nrow=nrow(X))%*%colMeans(X)
  output=X-m_mean
  return(output)}

#2. Covariance Matrix Funtion
cov_matrix<-function(X){
  m=ncol(X) #Number of columns of input matrix
  n=nrow(X) #Number of rows of input matrix
  
  #Compute X-E(X) for each element of the matrix
  m_mean=matrix(data=1, nrow=nrow(X))%*%colMeans(X)
  D=X-m_mean
  
  #Computing each element of the matrix
  C<-((n-1)^-1)*t(D)%*%D #Equivalent to computing [X-E(X)]^2
  return(C)
}

#3. Whitening Functions
###3.1 PCA Whitening
whiten_PCA<-function(X){
  #Centering
  X_c=centering_icaR(X)
  
  #Compute Covariance Matrix
  S=cov_matrix(X_c)
  
  #Eigen Value Decomposition
  D=eigen(S)$values
  for (i in 1:length(D)){
    if(D[i]<0){D[i]=0}
  }
  U=eigen(S)$vectors
  #Multiply by the sign of the eigen values
  U=sweep(U, 2, sign(diag(U)), "*") 
  
  #Regularization step
  epsilon = 1e-10 #To prevent NaNs
  
  #PCA Whitening Matrix
  W=(diag((D+epsilon)^(-1/2)))%*%t(U)
  
  #Xpca
  Xtilde=X%*%t(W)
  output=Xtilde
  return (output)
  }

##3.2 Whiten ZCA
whiten_ZCA<-function(X){
  #Centering
  X_c=centering_icaR(X)
  
  #Compute Covariance Matrix
  S=cov_matrix(X_c)
  
  #Eigen Value Decomposition
  D=eigen(S)$values
  for (i in 1:length(D)){
    if(D[i]<0){D[i]=0}
  }
  U=eigen(S)$vectors
  #Multiply by the sign of the eigen values
  U=sweep(U, 2, sign(diag(U)), "*") 
  
  #Regularization step
  epsilon = 1e-10 #To prevent NaNs
  
  #ZCA Whitening Matrix
  W=U%*%(diag((D+epsilon)^(-1/2)))%*%t(U)
  
  #Xpzca
  Xtilde=X%*%t(W)
  output=Xtilde
  return (output)
}

###Non-quadratic sub-gaussian or super-gaussian functions###
#fsup and its derivatives
fsup<-function(x,a){
  output=((1/a)*log(cosh(a*x)))-((a/2)*x^2)
  return(output)
}

fsup_1<-function(x,a){
  output=tanh(a*x)-(a*x)
  return(output)
}

fsup_2<-function(x,a){
  output=(a*((sech(a*x))^2)-1)
  return(output)
}

#fsub and its derivatives
fsub<-function(x,b){
  output=(b/4)*(x^4)
  return(output)
}

fsub_1<-function(x,b){
  output=b*x^3
  return(output)
}

fsub_2<-function(x,b){
  output=3*b*x^2
  return(output)
}

###Closeness Metrics###
#g(y,ref) = MSE
g.MSE<-function(y,ref){
  n=length(y)
  output=((y-ref)^2)
  return(output)
}
g.MSE_1<-function(y,ref){
  n=length(y)
  output=2*(y-ref)
  return(output)
}
g.MSE_2<-function(y,ref){
  n=length(y)
  output=2*(ones(n,1))
  return(output)
}
#g(y,ref) = cor
g.cor<-function(y,ref){
  n=length(y)
  output=(1/n)*((y*ref)^2)
  return(output)
}
g.cor_1<-function(y,ref){
  n=length(y)
  output=-2*(((1/n)*ref)/((1/n)*(y*ref)^3))
  return(output)
}
g.cor_2<-function(y,ref){
  n=length(y)
  output=6*((1/n)*(ref)^2)/((1/n)*(y*ref)^4)
  return(output)
}

######4. Lu & Rajapakse ICA with reference ICA-R Prototype######
ICA_R<-function(X, ref, threshold, learningRate, mu, lambda, gamma, a, b, ro, maxIter, OverValue, print = TRUE){
###i. Define dimensions of y (ICs)
ICnum=ncol(X)
IClen=nrow(X)
###ii.a Initialize parameters
mu=mu0
lambda0=lambda
gamma=1
a=a
b=b
ro=ro
ref=ref
v_gauss=rnorm(IClen, mean=0, sd=1)
###ii.a Initialize Weights
pInvX=round(ginv(X),14)
w0=(pInvX%*%ref)
w=w0/norm(w0, type=c("F"))
oldw=w
###ii.b Standarize Reference Signal
ref.stand=(ref-mean(ref))/sd(ref)
###iii. Center and Whiten data matrix X
X.center=centering_icaR(X)
X.whiten=whiten_PCA(X.center)
#a. Compute Moore-Penrose generalized inverse of X
pInvRxx=round(ginv(cov_matrix(X.whiten)),14)
#Define flow control parameters
flag=1
loop=1

while(flag==1){
###Compute initial y
y=X.whiten%*%w
y.stand=(y-mean(y))/sd(y)

  
###iv. Compute ro.hat
ro.hat=2*ro*(mean(fsup(a,y)-fsup(a,v_gauss)))

###v. Compute d using g'' as the correlation
d=mean(mu*g.MSE_2(y.stand,ref.stand))+(8*lambda)-mean(ro.hat*fsup_2(a,y))

###vi. Update mu/lambda multipliers
mu=max(-mu, gamma*(g.MSE(y.stand,ref.stand)-(threshold*(1-exp(-loop)))))
lambda=lambda+(gamma*(mean(y^2)-1)^2)

###vii. Compute phi and psi ###TO DO: Fix Expectation
phi=(1/IClen)*(t((-ro.hat*fsup_1(a,y)))%*%X.whiten)*(1/d)
psi=(1/IClen)*(t(mu*g.MSE_1(y.stand,ref.stand))%*%X.whiten)*(1/d)
omega=(4*lambda)*(mean(y^2)-1)*((1/IClen)*(t(y)%*%X.whiten))*(1/d)

###viii. Update/Normalize Weights
w=w-t((learningRate*(phi+psi+omega))%*%pInvRxx)
w=w/norm(w, type=c("F"))

###ix. Convergence criteria
wchange=1-as.numeric(abs(t(oldw)%*%w))
if(wchange<OverValue){
  flag=0
  print(sprintf("Converged after %d iteration",loop))
}

if(loop>=maxIter){
  flag=0
  print(sprintf("After %d iteration, still cannot converge",loop))

}
print(wchange)
oldw=w
loop=loop+1

}

y=X.whiten%*%w
print("End of cICA algorithm")

output=list(y,w)
return(output) 

}
#########################################################################################################
######5. Zhillin Zhang ICA with Reference ######
Zhang_cICA<-function(X, ref, threshold, learningRate, mu, lambda0, a, b, gamma, maxIter, OverValue, print=TRUE){
###i. Define dimensions of y (ICs)
ICnum=ncol(X)
IClen=nrow(X)
###ii.a Initialize parameters
mu=mu0
lambda=lambda0
gamma=1
rou=1
a=a
b=b
ref=ref
###ii.a Initialize Weights
w0=matrix(runif(ICnum))
w=w0/norm(w0, type=c("F"))
oldw=w
###ii.b Standarize Reference Signal
ref.stand=(ref-mean(ref))/sd(ref)
###iii. Center and Whiten data matrix X
X.center=centering_icaR(X)  
X.whiten=whiten_PCA(X.center)
#a. Compute Moore-Penrose generalized inverse of X
pInvRxx=round(ginv(cov_matrix(X.whiten)),14)
#Define flow control parameters
flag=1
loop=1
  
#While loop
while (flag==1) {
#1.Output using current iteration and standarize y
y=X.whiten%*%w
y.stand=(y-mean(y))/sd(y)
#2. Generate Gaussian Random Matrix and compute rou.hat
v_gauss=rnorm(IClen, mean=0, 1) 
rou.hat=rou*sign((mean(fsup(a,y)-fsup(a,v_gauss))))
    
#3. Gamma1 and Gamma2
Gamma1=((rou.hat*(t(fsup_1(a,y))%*%X.whiten))/IClen) - ((mu/2)*(t(g.MSE_1(y.stand, ref.stand))%*%X.whiten)/(IClen)) - ((lambda*t(y)%*%X.whiten)/(IClen))
Gamma2=((rou.hat*mean(fsup_2(a,y)))) - ((mu/2)*(mean(g.MSE_2(y.stand, ref.stand)))) - lambda

#4. Update Equations for mu and Lambda
mu=max(0, mu + (gamma * (g.MSE(y.stand,ref.stand)-(threshold*(1-exp(-loop))))))
lambda =lambda+(gamma*(mean(y.stand^2)-1)^2)

#5. Update/Normalize Weight Vector
w=w-(t(learningRate*(Gamma1*(1/Gamma2))%*%pInvRxx))
w=w/norm(w, type=c("F"))

#6.Check for Convergence
wchange=1-as.numeric(abs(t(oldw)%*%w))
if(wchange<OverValue){
  flag=0
  print(sprintf("Converged after %d iteration",loop))}

if(loop>=maxIter){
  flag=0
  print(sprintf("After %d iteration, still cannot converge",loop))
  
}
print(wchange)
oldw=w
loop=loop+1

}

y=X.whiten%*%w
print("End of cICA algorithm")

output=list(y,w)
return(output) 
}
#########################################################################################################

######6. Lin fastICA with Reference ######
Lin_fastICAR<-function(X, ref, threshold, learningRate, mu, a, b, ro, gamma, maxIter, OverValue, print=TRUE){
###i. Define dimensions of y (ICs)
ICnum=ncol(X)
IClen=nrow(X)
###ii.a Initialize parameters
mu=mu0
lambda=lambda0
gamma=gamma
rou=ro
a=a
b=b
ref=ref
###ii.a Initialize Weights
pInvX=round(ginv(X),14)
w0=(pInvX%*%ref)
w=w0/norm(w0, type=c("F"))
oldw=w
###ii.b Standarize Reference Signal
ref.stand=(ref-mean(ref))/sd(ref)
###iii. Center and Whiten data matrix X
X.center=centering_icaR(X)  
X.whiten=whiten_PCA(X.center)
#a. Compute Moore-Penrose generalized inverse of X
pInvRxx=round(ginv(cov_matrix(X.whiten)),14)
#Define flow control parameters
flag=1
loop=1
  
#While loop
while (flag==1) {
    
#1.Output using current iteration and standarize y
y=X.whiten%*%w
y.stand=(y-mean(y))/sd(y)
#2. Generate Gaussian Random Matrix and compute rou.hat
v_gauss=rnorm(IClen, mean=0, 1) 
rou.hat=(mean(fsup(a,y)-fsup(a,v_gauss)))
    
#3. L and delta
L=rou.hat*((1/IClen)*(t(fsup_1(a,y))%*%X.whiten)) - ((mu/2)*(t(g.MSE_1(y.stand, ref.stand))%*%X.whiten)/(IClen))
delta=((rou.hat*mean(fsup_2(a,y)))) - ((mu/2)*(mean(g.MSE_2(y.stand, ref.stand))))

#4. Update Equations for mu 
mu=max(0, mu + (gamma * (g.MSE(y.stand,ref.stand)-(threshold*(1-exp(-loop))))))

#5. Update/Normalize Weight Vector
w=w-(t(learningRate*L*(1/delta)))
w=w/norm(w, type=c("F"))
    
#6.Check for Convergence
wchange=1-as.numeric(abs(t(oldw)%*%w))
  if(wchange<OverValue){
    flag=0
    print(sprintf("Converged after %d iteration",loop))
    }
    
  if(loop>=maxIter){
    flag=0
    print(sprintf("After %d iteration, still cannot converge",loop))
  }

print(wchange)
oldw=w
loop=loop+1

}
  
y=X.whiten%*%w
print("End of cICA algorithm")
  
output=list(y,w)
return(output) 
}
#########################################################################################################



####7. Signal To Noise Ratio###
SNR<-function(s,c){
  Var.S=var(s)
  MSE=(sum(g.MSE(c,s)))/length(S[,1])
  SNR=10*(log10(Var.S/MSE))
  return(SNR)
}
################################

####8. Individual Performance Index ####
IPI<-function(w,A){
  #i. Compute Permutation Matrix
  Perm.Vec=t(w)%*%A
  
  #ii. Define dimensions and maximum permutation vector
  M=length(Perm.Vec)
  p.k=max(abs(Perm.Vec))
  
  #iii. Compute IPI
  s.IPI=0
  for (j in 1:M){
    p.j=abs(Perm.Vec[,j])
    aux=p.j/p.k
    s.IPI=s.IPI+aux
    #print(aux)
    #print(s.IPI)
  }
  IPI=s.IPI-1
  #i.v Return IPI
  return(IPI)
}
##########################################
###TEST 1###
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
S[,4] = rnorm(N)
S[,5] = rnorm(N)

A = matrix(runif(25), 5,5)
X = S%*%A

ref1=S[,1]
ref2 = cos(2*pi*f2*ts*k)

mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1


c11=ICA_R(X, ref1, threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
c12=ICA_R(X, ref2, threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
c2=Zhang_cICA(X, ref1, threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
c31=Lin_fastICAR(X, ref1, threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
c32=Lin_fastICAR(X, ref2, threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)

plot(S[,1], type="l")
plot(S[,2], type="l")
plot(c11[[1]], type="l")
plot(c12[[1]], type="l")
plot(c2[[1]], type="l")
plot(c31[[1]], type="l")
plot(c32[[1]], type="l")

SNR(S[,1], c11[[1]])
SNR(S[,1], c2[[1]])
SNR(S[,1], c31[[1]])


IPI(c11[[2]], A)
IPI(c2[[2]], A)
IPI(c31[[2]], A)


#####################################TEST 2#####################################T
#Test matrixes Simulation 2: Random Frequencies#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.kPCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()

#SIMULATION#

for (j in 1:1000){
#SIMULATION PARAMETERS#  
ts = 1e-4 # sampling period
f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
f2 = runif(N, min=0.001, max=0.1)/ts    
f3 = runif(N, min=0.001, max=0.1)/ts
S=matrix(0, N, 5)
  
S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k))
S[,2] = cos(2*pi*f2*ts*k)
S[,3] = cos(2*pi*f3*ts*k + 2)
S[,4] = rnorm(N)
S[,5] = rnorm(N)
  
A = matrix(runif(25), 5,5)
X = S%*%A
L=3 #Three relevant signals
#cICA parameters
ref1=S[,1]
ref2=S[,2]
ref3=S[,3]
ref=data.frame(ref1, ref2, ref3)
#Run fastICA
ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
#Run Kernel PCA using Radial-Basis Function
kpc=kpca(scale(X), kernel="rbfdot", features=3)

for (i in 1:L){
  #Run cICA Algorithms
  c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
  #Compute SNR
  SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
  SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
  SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
  SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
  SNR.df[i,5]=as.numeric(SNR(pcv(kpc)[,i], S[,i]))
  
  #Compute IPI (Ad-hoc measurement for ICA framework)
  IPI.df[i,1]=IPI(c1[[2]], A)
  IPI.df[i,2]=IPI(c2[[2]], A)
  IPI.df[i,3]=IPI(c3[[2]], A)
  IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)

}
colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "kPCA")
colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")

#Store Simulation Results
#SNR
aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
aux.SNR.kPCA=rbind.data.frame(aux.SNR.kPCA, SNR.df[,5])
#IPI
aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.kPCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#########SNR RESULTS#######
for (i in 1:dim(aux.SNR.c1)[2]){
  #c1
  phistSNR1<-ggplot(aux.SNR.c1, aes(x=aux.SNR.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_LR_", colnames(aux.SNR.c1)[i],
                    ".PNG", sep="")
  ggsave(phistSNR1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistSNR2<-ggplot(aux.SNR.c2, aes(x=aux.SNR.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_Zhang_", colnames(aux.SNR.c2)[i],
                    ".PNG", sep="")
  ggsave(phistSNR2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistSNR3<-ggplot(aux.SNR.c3, aes(x=aux.SNR.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_Lin_", colnames(aux.SNR.c3)[i],
                    ".PNG", sep="")
  ggsave(phistSNR3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistSNR4<-ggplot(aux.SNR.fastICA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_fastICA_", colnames(aux.SNR.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR4, file=file_name, width = 14, height = 10, units = "cm")
  #kPCA
  phistSNR5<-ggplot(aux.SNR.kPCA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_kPCA_", colnames(aux.SNR.kPCA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR5, file=file_name, width = 14, height = 10, units = "cm")
}

#SNR side-to-side comparison
for (i in 1:3){
a1=data.frame("LR",aux.SNR.c1[,i])
colnames(a1)=c("Algorithm", "SNR")
a2=data.frame("Zhang",aux.SNR.c2[,i])
colnames(a2)=c("Algorithm", "SNR")
a3=data.frame("Lin",aux.SNR.c3[,i])
colnames(a3)=c("Algorithm", "SNR")
a4=data.frame("fastICA",aux.SNR.fastICA[,i])
colnames(a4)=c("Algorithm", "SNR")
a5=data.frame("kPCA",aux.SNR.kPCA[,i])
colnames(a5)=c("Algorithm", "SNR")
a=rbind(a1, a2, a3, a4, a5)
p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
  xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/SNR/SNR_benchmark_", 
                  colnames(aux.SNR.c1)[i],".PNG", sep="")
ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute SNR Median
m.SNR.c1=lapply(aux.SNR.c1, median)
m.SNR.c2=lapply(aux.SNR.c2, median)
m.SNR.c3=lapply(aux.SNR.c3, median)
m.SNR.fastICA=lapply(aux.SNR.fastICA, median)
m.SNR.kPCA=lapply(aux.SNR.kPCA, median)

#Compute SNR Quantiles
q.SNR.c1=lapply(aux.SNR.c1, quantile,  probs=c(0.025,0.975))
q.SNR.c2=lapply(aux.SNR.c2, quantile,  probs=c(0.025,0.975))
q.SNR.c3=lapply(aux.SNR.c3, quantile,  probs=c(0.025,0.975))
q.SNR.fastICA=lapply(aux.SNR.fastICA, quantile,  probs=c(0.025,0.975))
q.SNR.kPCA=lapply(aux.SNR.kPCA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
kPCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.kPCA[,i]
  bootCI.SNR.kPCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=10000))
  kPCAResults=rbind.data.frame(kPCAResults,data.frame(bootCI.SNR.kPCA$basic[4], median(x5), bootCI.SNR.kPCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(kPCAResults)=c("kPCA.lowerCI", "kPCA.Median.SNR", "kPCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, kPCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/SNR_Results_Sim1.csv", row.names = FALSE)

##################################################################################
#########IPI RESULTS#######
for (i in 1:dim(aux.IPI.c1)[2]){
  #c1
  phistIPI1<-ggplot(aux.IPI.c1, aes(x=aux.IPI.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/IPI/IPI_LR_", colnames(aux.IPI.c1)[i],
                    ".PNG", sep="")
  ggsave(phistIPI1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistIPI2<-ggplot(aux.IPI.c2, aes(x=aux.IPI.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/IPI/IPI_Zhang_", colnames(aux.IPI.c2)[i],
                    ".PNG", sep="")
  ggsave(phistIPI2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistIPI3<-ggplot(aux.IPI.c3, aes(x=aux.IPI.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/IPI/IPI_Lin_", colnames(aux.IPI.c3)[i],
                    ".PNG", sep="")
  ggsave(phistIPI3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistIPI4<-ggplot(aux.IPI.fastICA, aes(x=aux.IPI.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/IPI/IPI_fastICA_", colnames(aux.IPI.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistIPI4, file=file_name, width = 14, height = 10, units = "cm")
}

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 1 Noiseless/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute IPI Median
m.IPI.c1=lapply(aux.IPI.c1, median)
m.IPI.c2=lapply(aux.IPI.c2, median)
m.IPI.c3=lapply(aux.IPI.c3, median)
m.IPI.fastICA=lapply(aux.IPI.fastICA, median)
#Compute IPI Quantiles
q.IPI.c1=lapply(aux.IPI.c1, quantile,  probs=c(0.025,0.975))
q.IPI.c2=lapply(aux.IPI.c2, quantile,  probs=c(0.025,0.975))
q.IPI.c3=lapply(aux.IPI.c3, quantile,  probs=c(0.025,0.975))
q.IPI.fastICA=lapply(aux.IPI.fastICA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/IPI_Results_Sim1.csv", row.names = FALSE)


#####################################TEST 3#####################################T
#Test matrixes Simulation 3: Random Frequencies and White Noise#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.kPCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()

#SIMULATION#

for (j in 1:1000){
  #SIMULATION PARAMETERS#  
  ts = 1e-4 # sampling period
  f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
  f2 = runif(N, min=0.001, max=0.1)/ts    
  f3 = runif(N, min=0.001, max=0.1)/ts
  S=matrix(0, N, 5)
  
  S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k)) 
  S[,2] = cos(2*pi*f2*ts*k)
  S[,3] = cos(2*pi*f3*ts*k + 2)
  S[,4] = rnorm(N)
  S[,5] = rnorm(N)
  
  A = matrix(runif(25), 5,5)
  X = S%*%A
  L=3 #Three relevant signals
  #cICA parameters
  ref1=S[,1] + rnorm(N, 0, 1)
  ref2=S[,2] + rnorm(N, 0, 1)
  ref3=S[,3] + rnorm(N, 0, 1)
  ref=data.frame(ref1, ref2, ref3)
  #Run fastICA
  ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
  #Run Kernel PCA using Radial-Basis Function
  kpc=kpca(scale(X), kernel="rbfdot", features=3)
  
  for (i in 1:L){
    #Run cICA Algorithms
    c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    #Compute SNR
    SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
    SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
    SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
    SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
    SNR.df[i,5]=as.numeric(SNR(pcv(kpc)[,i], S[,i]))
    
    #Compute IPI (Ad-hoc measurement for ICA framework)
    IPI.df[i,1]=IPI(c1[[2]], A)
    IPI.df[i,2]=IPI(c2[[2]], A)
    IPI.df[i,3]=IPI(c3[[2]], A)
    IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)
    
  }
  colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "kPCA")
  colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
  
  #Store Simulation Results
  #SNR
  aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
  aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
  aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
  aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
  aux.SNR.kPCA=rbind.data.frame(aux.SNR.kPCA, SNR.df[,5])
  #IPI
  aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
  aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
  aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
  aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.kPCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#########SNR RESULTS#######
for (i in 1:dim(aux.SNR.c1)[2]){
  #c1
  phistSNR1<-ggplot(aux.SNR.c1, aes(x=aux.SNR.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/SNR/SNR_LR_", colnames(aux.SNR.c1)[i],
                    ".PNG", sep="")
  ggsave(phistSNR1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistSNR2<-ggplot(aux.SNR.c2, aes(x=aux.SNR.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/SNR/SNR_Zhang_", colnames(aux.SNR.c2)[i],
                    ".PNG", sep="")
  ggsave(phistSNR2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistSNR3<-ggplot(aux.SNR.c3, aes(x=aux.SNR.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/SNR/SNR_Lin_", colnames(aux.SNR.c3)[i],
                    ".PNG", sep="")
  ggsave(phistSNR3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistSNR4<-ggplot(aux.SNR.fastICA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/SNR/SNR_fastICA_", colnames(aux.SNR.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR4, file=file_name, width = 14, height = 10, units = "cm")
  #kPCA
  phistSNR5<-ggplot(aux.SNR.kPCA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/SNR/SNR_kPCA_", colnames(aux.SNR.kPCA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR5, file=file_name, width = 14, height = 10, units = "cm")
}

#SNR side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.SNR.c1[,i])
  colnames(a1)=c("Algorithm", "SNR")
  a2=data.frame("Zhang",aux.SNR.c2[,i])
  colnames(a2)=c("Algorithm", "SNR")
  a3=data.frame("Lin",aux.SNR.c3[,i])
  colnames(a3)=c("Algorithm", "SNR")
  a4=data.frame("fastICA",aux.SNR.fastICA[,i])
  colnames(a4)=c("Algorithm", "SNR")
  a5=data.frame("kPCA",aux.SNR.kPCA[,i])
  colnames(a5)=c("Algorithm", "SNR")
  a=rbind(a1, a2, a3, a4, a5)
  p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/SNR/SNR_benchmark_", 
                    colnames(aux.SNR.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute SNR Median
m.SNR.c1=lapply(aux.SNR.c1, median)
m.SNR.c2=lapply(aux.SNR.c2, median)
m.SNR.c3=lapply(aux.SNR.c3, median)
m.SNR.fastICA=lapply(aux.SNR.fastICA, median)
m.SNR.kPCA=lapply(aux.SNR.kPCA, median)

#Compute SNR Quantiles
q.SNR.c1=lapply(aux.SNR.c1, quantile,  probs=c(0.025,0.975))
q.SNR.c2=lapply(aux.SNR.c2, quantile,  probs=c(0.025,0.975))
q.SNR.c3=lapply(aux.SNR.c3, quantile,  probs=c(0.025,0.975))
q.SNR.fastICA=lapply(aux.SNR.fastICA, quantile,  probs=c(0.025,0.975))
q.SNR.kPCA=lapply(aux.SNR.kPCA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
kPCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.kPCA[,i]
  bootCI.SNR.kPCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=10000))
  kPCAResults=rbind.data.frame(kPCAResults,data.frame(bootCI.SNR.kPCA$basic[4], median(x5), bootCI.SNR.kPCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(kPCAResults)=c("kPCA.lowerCI", "kPCA.Median.SNR", "kPCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, kPCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/SNR_Results_Sim2.csv", row.names = FALSE)

##################################################################################
#########IPI RESULTS#######
for (i in 1:dim(aux.IPI.c1)[2]){
  #c1
  phistIPI1<-ggplot(aux.IPI.c1, aes(x=aux.IPI.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/IPI/IPI_LR_", colnames(aux.IPI.c1)[i],
                    ".PNG", sep="")
  ggsave(phistIPI1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistIPI2<-ggplot(aux.IPI.c2, aes(x=aux.IPI.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/IPI/IPI_Zhang_", colnames(aux.IPI.c2)[i],
                    ".PNG", sep="")
  ggsave(phistIPI2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistIPI3<-ggplot(aux.IPI.c3, aes(x=aux.IPI.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/IPI/IPI_Lin_", colnames(aux.IPI.c3)[i],
                    ".PNG", sep="")
  ggsave(phistIPI3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistIPI4<-ggplot(aux.IPI.fastICA, aes(x=aux.IPI.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/IPI/IPI_fastICA_", colnames(aux.IPI.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistIPI4, file=file_name, width = 14, height = 10, units = "cm")
}

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 2 Gaussian Noise/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute IPI Median
m.IPI.c1=lapply(aux.IPI.c1, median)
m.IPI.c2=lapply(aux.IPI.c2, median)
m.IPI.c3=lapply(aux.IPI.c3, median)
m.IPI.fastICA=lapply(aux.IPI.fastICA, median)
#Compute IPI Quantiles
q.IPI.c1=lapply(aux.IPI.c1, quantile,  probs=c(0.025,0.975))
q.IPI.c2=lapply(aux.IPI.c2, quantile,  probs=c(0.025,0.975))
q.IPI.c3=lapply(aux.IPI.c3, quantile,  probs=c(0.025,0.975))
q.IPI.fastICA=lapply(aux.IPI.fastICA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/IPI_Results_Sim2.csv", row.names = FALSE)


#####################################TEST 4#####################################T
#Test matrixes Simulation 4: Random Frequencies and Gaussian Noise with same mean and variance as the source signals#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.kPCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()

#SIMULATION#

for (j in 1:1000){
  #SIMULATION PARAMETERS#  
  ts = 1e-4 # sampling period
  f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
  f2 = runif(N, min=0.001, max=0.1)/ts    
  f3 = runif(N, min=0.001, max=0.1)/ts
  S=matrix(0, N, 5)
  
  S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k)) 
  S[,2] = cos(2*pi*f2*ts*k)
  S[,3] = cos(2*pi*f3*ts*k + 2)
  S[,4] = rnorm(N)
  S[,5] = rnorm(N)
  
  A = matrix(runif(25), 5,5)
  X = S%*%A
  L=3 #Three relevant signals
  #cICA parameters
  ref1=S[,1] + rnorm(N, mean(S[,1]), sd(S[,1]))
  ref2=S[,2] + rnorm(N, mean(S[,2]), sd(S[,2]))
  ref3=S[,3] + rnorm(N, mean(S[,3]), sd(S[,3]))
  ref=data.frame(ref1, ref2, ref3)
  #Run fastICA
  ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
  #Run Kernel PCA using Radial-Basis Function
  kpc=kpca(scale(X), kernel="rbfdot", features=3)
  
  for (i in 1:L){
    #Run cICA Algorithms
    c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    #Compute SNR
    SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
    SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
    SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
    SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
    SNR.df[i,5]=as.numeric(SNR(pcv(kpc)[,i], S[,i]))
    
    #Compute IPI (Ad-hoc measurement for ICA framework)
    IPI.df[i,1]=IPI(c1[[2]], A)
    IPI.df[i,2]=IPI(c2[[2]], A)
    IPI.df[i,3]=IPI(c3[[2]], A)
    IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)
    
  }
  colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "kPCA")
  colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
  
  #Store Simulation Results
  #SNR
  aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
  aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
  aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
  aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
  aux.SNR.kPCA=rbind.data.frame(aux.SNR.kPCA, SNR.df[,5])
  #IPI
  aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
  aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
  aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
  aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.kPCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#########SNR RESULTS#######
for (i in 1:dim(aux.SNR.c1)[2]){
  #c1
  phistSNR1<-ggplot(aux.SNR.c1, aes(x=aux.SNR.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_LR_", colnames(aux.SNR.c1)[i],
                    ".PNG", sep="")
  ggsave(phistSNR1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistSNR2<-ggplot(aux.SNR.c2, aes(x=aux.SNR.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_Zhang_", colnames(aux.SNR.c2)[i],
                    ".PNG", sep="")
  ggsave(phistSNR2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistSNR3<-ggplot(aux.SNR.c3, aes(x=aux.SNR.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_Lin_", colnames(aux.SNR.c3)[i],
                    ".PNG", sep="")
  ggsave(phistSNR3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistSNR4<-ggplot(aux.SNR.fastICA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_fastICA_", colnames(aux.SNR.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR4, file=file_name, width = 14, height = 10, units = "cm")
  #kPCA
  phistSNR5<-ggplot(aux.SNR.kPCA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_kPCA_", colnames(aux.SNR.kPCA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR5, file=file_name, width = 14, height = 10, units = "cm")
}

#SNR side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.SNR.c1[,i])
  colnames(a1)=c("Algorithm", "SNR")
  a2=data.frame("Zhang",aux.SNR.c2[,i])
  colnames(a2)=c("Algorithm", "SNR")
  a3=data.frame("Lin",aux.SNR.c3[,i])
  colnames(a3)=c("Algorithm", "SNR")
  a4=data.frame("fastICA",aux.SNR.fastICA[,i])
  colnames(a4)=c("Algorithm", "SNR")
  a5=data.frame("kPCA",aux.SNR.kPCA[,i])
  colnames(a5)=c("Algorithm", "SNR")
  a=rbind(a1, a2, a3, a4, a5)
  p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/SNR/SNR_benchmark_", 
                    colnames(aux.SNR.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute SNR Median
m.SNR.c1=lapply(aux.SNR.c1, median)
m.SNR.c2=lapply(aux.SNR.c2, median)
m.SNR.c3=lapply(aux.SNR.c3, median)
m.SNR.fastICA=lapply(aux.SNR.fastICA, median)
m.SNR.kPCA=lapply(aux.SNR.kPCA, median)

#Compute SNR Quantiles
q.SNR.c1=lapply(aux.SNR.c1, quantile,  probs=c(0.025,0.975))
q.SNR.c2=lapply(aux.SNR.c2, quantile,  probs=c(0.025,0.975))
q.SNR.c3=lapply(aux.SNR.c3, quantile,  probs=c(0.025,0.975))
q.SNR.fastICA=lapply(aux.SNR.fastICA, quantile,  probs=c(0.025,0.975))
q.SNR.kPCA=lapply(aux.SNR.kPCA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
kPCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.kPCA[,i]
  bootCI.SNR.kPCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=10000))
  kPCAResults=rbind.data.frame(kPCAResults,data.frame(bootCI.SNR.kPCA$basic[4], median(x5), bootCI.SNR.kPCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(kPCAResults)=c("kPCA.lowerCI", "kPCA.Median.SNR", "kPCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, kPCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/SNR_Results_Sim3.csv", row.names = FALSE)

##################################################################################
#########IPI RESULTS#######
for (i in 1:dim(aux.IPI.c1)[2]){
  #c1
  phistIPI1<-ggplot(aux.IPI.c1, aes(x=aux.IPI.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/IPI/IPI_LR_", colnames(aux.IPI.c1)[i],
                    ".PNG", sep="")
  ggsave(phistIPI1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistIPI2<-ggplot(aux.IPI.c2, aes(x=aux.IPI.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/IPI/IPI_Zhang_", colnames(aux.IPI.c2)[i],
                    ".PNG", sep="")
  ggsave(phistIPI2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistIPI3<-ggplot(aux.IPI.c3, aes(x=aux.IPI.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/IPI/IPI_Lin_", colnames(aux.IPI.c3)[i],
                    ".PNG", sep="")
  ggsave(phistIPI3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistIPI4<-ggplot(aux.IPI.fastICA, aes(x=aux.IPI.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/IPI/IPI_fastICA_", colnames(aux.IPI.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistIPI4, file=file_name, width = 14, height = 10, units = "cm")
}

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 3 Gaussian Noise/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute IPI Median
m.IPI.c1=lapply(aux.IPI.c1, median)
m.IPI.c2=lapply(aux.IPI.c2, median)
m.IPI.c3=lapply(aux.IPI.c3, median)
m.IPI.fastICA=lapply(aux.IPI.fastICA, median)
#Compute IPI Quantiles
q.IPI.c1=lapply(aux.IPI.c1, quantile,  probs=c(0.025,0.975))
q.IPI.c2=lapply(aux.IPI.c2, quantile,  probs=c(0.025,0.975))
q.IPI.c3=lapply(aux.IPI.c3, quantile,  probs=c(0.025,0.975))
q.IPI.fastICA=lapply(aux.IPI.fastICA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/IPI_Results_Sim3.csv", row.names = FALSE)




#####################################TEST 5#####################################
library("aTSA", lib.loc="~/R/win-library/3.5")
library("tseries", lib.loc="~/R/win-library/3.5")
library("TSA", lib.loc="~/R/win-library/3.5")
library("forecast", lib.loc="~/R/win-library/3.5")

#Test matrixes Simulation 5: Random Frequencies and Seasonally Adjusted Reference Signals#
library("fastICA", lib.loc="~/R/win-library/3.5")
library("kernlab", lib.loc="~/R/win-library/3.5")
#SIMULATION HYPERPARAMETERS
N=300
k=seq(1:N)
ts = 1e-4 # sampling period
#cICA parameters
mu0 = 1
lambda0 = 1
gamma0 = 1
learningRate0 = 1
OverValue0=1e-6  
maxIter0 = 200000
threshold0 = 1.5
a0=1
b0=1
ro0=1
#Empty Data Frames
aux.SNR.c1=data.frame()
aux.SNR.c2=data.frame()
aux.SNR.c3=data.frame()
aux.SNR.fastICA=data.frame()
aux.SNR.kPCA=data.frame()
SNR.df=data.frame()
IPI.df=data.frame()
aux.IPI.c1=data.frame()
aux.IPI.c2=data.frame()
aux.IPI.c3=data.frame()
aux.IPI.fastICA=data.frame()


#SIMULATION#

for (j in 1:1000){
  #SIMULATION PARAMETERS#  
  ts = 1e-4 # sampling period
  f1 = runif(N, min=0.001, max=0.1)/ts    # true frequency
  f2 = runif(N, min=0.001, max=0.1)/ts    
  f3 = runif(N, min=0.001, max=0.1)/ts
  S=matrix(0, N, 5)
  
  S[,1] = sin(2*pi*f1*ts*k +6*cos(2*pi*200*ts*k)) 
  S[,2] = cos(2*pi*f2*ts*k)
  S[,3] = cos(2*pi*f3*ts*k + 2)
  S[,4] = rnorm(N)
  S[,5] = rnorm(N)
  
  A = matrix(runif(25), 5,5)
  X = S%*%A
  L=3 #Three relevant signals
  
  #cICA parameters
  p1=spec.pgram(X[,1], plot=FALSE, taper=0.1, demean=TRUE)
  freq1=p1$freq[which(p1$spec>=sort(p1$spec, decreasing = TRUE)[3], arr.ind = TRUE)]
  ref1=sin(freq1[1]*k)
  p2=spec.pgram(X[,2], plot=FALSE, taper=0.1, demean=TRUE)
  freq2=p2$freq[which(p2$spec>=sort(p2$spec, decreasing = TRUE)[3], arr.ind = TRUE)]
  ref2=cos(freq2[1]*k)
  p3=spec.pgram(X[,3], plot=FALSE, taper=0.1, demean=TRUE)
  freq3=p3$freq[which(p3$spec>=sort(p3$spec, decreasing = TRUE)[3], arr.ind = TRUE)]
  ref3=cos(freq3[1]*k)
  ref=data.frame(ref1, ref2, ref3)
  #Run fastICA
  ic=fastICA(X, L, alg.typ = "deflation", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = maxIter0, tol = 0.0001, verbose = FALSE)
  #Run Kernel PCA using Radial-Basis Function
  kpc=kpca(scale(X), kernel="rbfdot", features=3)
  
  for (i in 1:L){
    #Run cICA Algorithms
    c1=ICA_R(X, ref[,i], threshold0, learningRate0, mu0, lambda0,  gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c2=Zhang_cICA(X, ref[,i], threshold0, learningRate0, mu0, lambda0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    c3=Lin_fastICAR(X, ref[,i], threshold0, learningRate0, mu0, gamma0, a0, b0, ro0, maxIter0, OverValue0)
    #Compute SNR
    SNR.df[i,1]=as.numeric(SNR(c1[[1]], S[,i]))
    SNR.df[i,2]=as.numeric(SNR(c2[[1]], S[,i]))
    SNR.df[i,3]=as.numeric(SNR(c3[[1]], S[,i]))
    SNR.df[i,4]=as.numeric(SNR(ic$S[,i],S[,i]))
    SNR.df[i,5]=as.numeric(SNR(pcv(kpc)[,i], S[,i]))
    
    #Compute IPI (Ad-hoc measurement for ICA framework)
    IPI.df[i,1]=IPI(c1[[2]], A)
    IPI.df[i,2]=IPI(c2[[2]], A)
    IPI.df[i,3]=IPI(c3[[2]], A)
    IPI.df[i,4]=IPI(as.matrix(ic$A[i,]),A)
    
  }
  colnames(SNR.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA", "kPCA")
  colnames(IPI.df)=c("LR.cICA", "Zhang.cICA", "Lin.cICA", "fastICA")
  
  #Store Simulation Results
  #SNR
  aux.SNR.c1=rbind.data.frame(aux.SNR.c1, SNR.df[,1])
  aux.SNR.c2=rbind.data.frame(aux.SNR.c2, SNR.df[,2])
  aux.SNR.c3=rbind.data.frame(aux.SNR.c3, SNR.df[,3])
  aux.SNR.fastICA=rbind.data.frame(aux.SNR.fastICA, SNR.df[,4])
  aux.SNR.kPCA=rbind.data.frame(aux.SNR.kPCA, SNR.df[,5])
  #IPI
  aux.IPI.c1=rbind.data.frame(aux.IPI.c1, IPI.df[,1])
  aux.IPI.c2=rbind.data.frame(aux.IPI.c2, IPI.df[,2])
  aux.IPI.c3=rbind.data.frame(aux.IPI.c3, IPI.df[,3])
  aux.IPI.fastICA=rbind.data.frame(aux.IPI.fastICA, IPI.df[,4])
}

#Rename dataframes
#SNR
colnames(aux.SNR.c1)=c("S1", "S2", "S3")
colnames(aux.SNR.c2)=c("S1", "S2", "S3")
colnames(aux.SNR.c3)=c("S1", "S2", "S3")
colnames(aux.SNR.fastICA)=c("S1", "S2", "S3")
colnames(aux.SNR.kPCA)=c("S1", "S2", "S3")
#IPI
colnames(aux.IPI.c1)=c("S1", "S2", "S3")
colnames(aux.IPI.c2)=c("S1", "S2", "S3")
colnames(aux.IPI.c3)=c("S1", "S2", "S3")
colnames(aux.IPI.fastICA)=c("S1", "S2", "S3")


####Plotting Results####
library("boot", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
#########SNR RESULTS#######
for (i in 1:dim(aux.SNR.c1)[2]){
  #c1
  phistSNR1<-ggplot(aux.SNR.c1, aes(x=aux.SNR.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_LR_", colnames(aux.SNR.c1)[i],
                    ".PNG", sep="")
  ggsave(phistSNR1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistSNR2<-ggplot(aux.SNR.c2, aes(x=aux.SNR.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_Zhang_", colnames(aux.SNR.c2)[i],
                    ".PNG", sep="")
  ggsave(phistSNR2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistSNR3<-ggplot(aux.SNR.c3, aes(x=aux.SNR.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_Lin_", colnames(aux.SNR.c3)[i],
                    ".PNG", sep="")
  ggsave(phistSNR3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistSNR4<-ggplot(aux.SNR.fastICA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_fastICA_", colnames(aux.SNR.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR4, file=file_name, width = 14, height = 10, units = "cm")
  #kPCA
  phistSNR5<-ggplot(aux.SNR.kPCA, aes(x=aux.SNR.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_kPCA_", colnames(aux.SNR.kPCA)[i],
                    ".PNG", sep="")
  ggsave(phistSNR5, file=file_name, width = 14, height = 10, units = "cm")
}

#SNR side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.SNR.c1[,i])
  colnames(a1)=c("Algorithm", "SNR")
  a2=data.frame("Zhang",aux.SNR.c2[,i])
  colnames(a2)=c("Algorithm", "SNR")
  a3=data.frame("Lin",aux.SNR.c3[,i])
  colnames(a3)=c("Algorithm", "SNR")
  a4=data.frame("fastICA",aux.SNR.fastICA[,i])
  colnames(a4)=c("Algorithm", "SNR")
  a5=data.frame("kPCA",aux.SNR.kPCA[,i])
  colnames(a5)=c("Algorithm", "SNR")
  a=rbind(a1, a2, a3, a4, a5)
  p=ggplot(a, aes(x=SNR, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Signal-to-Noise Ratio") + ylab("Frequency Count")+labs(title = colnames(aux.SNR.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/SNR/SNR_benchmark_", 
                    colnames(aux.SNR.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute SNR Median
m.SNR.c1=lapply(aux.SNR.c1, median)
m.SNR.c2=lapply(aux.SNR.c2, median)
m.SNR.c3=lapply(aux.SNR.c3, median)
m.SNR.fastICA=lapply(aux.SNR.fastICA, median)
m.SNR.kPCA=lapply(aux.SNR.kPCA, median)

#Compute SNR Quantiles
q.SNR.c1=lapply(aux.SNR.c1, quantile,  probs=c(0.025,0.975))
q.SNR.c2=lapply(aux.SNR.c2, quantile,  probs=c(0.025,0.975))
q.SNR.c3=lapply(aux.SNR.c3, quantile,  probs=c(0.025,0.975))
q.SNR.fastICA=lapply(aux.SNR.fastICA, quantile,  probs=c(0.025,0.975))
q.SNR.kPCA=lapply(aux.SNR.kPCA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
kPCAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.SNR.c1[,i]
  bootCI.SNR.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.SNR.c1$basic[4], median(x1), bootCI.SNR.c1$basic[5]))
  #Zhang
  x2=aux.SNR.c2[,i]
  bootCI.SNR.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.SNR.c2$basic[4], median(x2), bootCI.SNR.c2$basic[5]))
  #Lin
  x3=aux.SNR.c3[,i]
  bootCI.SNR.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.SNR.c3$basic[4], median(x3), bootCI.SNR.c3$basic[5]))
  #fastICA
  x4=aux.SNR.fastICA[,i]
  bootCI.SNR.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.SNR.fastICA$basic[4], median(x4), bootCI.SNR.fastICA$basic[5]))
  #kPCA
  x5=aux.SNR.kPCA[,i]
  bootCI.SNR.kPCA=boot.ci(boot(x5, function(x,j) median(x[j]), R=10000))
  kPCAResults=rbind.data.frame(kPCAResults,data.frame(bootCI.SNR.kPCA$basic[4], median(x5), bootCI.SNR.kPCA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.SNR", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.SNR", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.SNR", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.SNR", "fastICA.upperCI")
colnames(kPCAResults)=c("kPCA.lowerCI", "kPCA.Median.SNR", "kPCA.upperCI")

SNR.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults, kPCAResults)
write.csv(SNR.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/SNR_Results_Sim4.csv", row.names = FALSE)

##################################################################################
#########IPI RESULTS#######
for (i in 1:dim(aux.IPI.c1)[2]){
  #c1
  phistIPI1<-ggplot(aux.IPI.c1, aes(x=aux.IPI.c1[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/IPI/IPI_LR_", colnames(aux.IPI.c1)[i],
                    ".PNG", sep="")
  ggsave(phistIPI1, file=file_name, width = 14, height = 10, units = "cm")
  #c2
  phistIPI2<-ggplot(aux.IPI.c2, aes(x=aux.IPI.c2[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/IPI/IPI_Zhang_", colnames(aux.IPI.c2)[i],
                    ".PNG", sep="")
  ggsave(phistIPI2, file=file_name, width = 14, height = 10, units = "cm")
  #c3
  phistIPI3<-ggplot(aux.IPI.c3, aes(x=aux.IPI.c3[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/IPI/IPI_Lin_", colnames(aux.IPI.c3)[i],
                    ".PNG", sep="")
  ggsave(phistIPI3, file=file_name, width = 14, height = 10, units = "cm")
  #fastICA
  phistIPI4<-ggplot(aux.IPI.fastICA, aes(x=aux.IPI.fastICA[,i])) + geom_histogram(aes(y=..density..), colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666") 
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/IPI/IPI_fastICA_", colnames(aux.IPI.fastICA)[i],
                    ".PNG", sep="")
  ggsave(phistIPI4, file=file_name, width = 14, height = 10, units = "cm")
}

#IPI side-to-side comparison
for (i in 1:3){
  a1=data.frame("LR",aux.IPI.c1[,i])
  colnames(a1)=c("Algorithm", "IPI")
  a2=data.frame("Zhang",aux.IPI.c2[,i])
  colnames(a2)=c("Algorithm", "IPI")
  a3=data.frame("Lin",aux.IPI.c3[,i])
  colnames(a3)=c("Algorithm", "IPI")
  a4=data.frame("fastICA",aux.IPI.fastICA[,i])
  colnames(a4)=c("Algorithm", "IPI")
  a=rbind(a1, a2, a3, a4)
  p=ggplot(a, aes(x=IPI, fill=Algorithm)) + geom_histogram(binwidth=.5, position="dodge") + 
    xlab("Individual Performance Index") + ylab("Frequency Count")+labs(title = colnames(aux.IPI.c1)[i])
  file_name = paste("C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results Plots/Sim 4 Rough Ref Signal/IPI/IPI_benchmark_", 
                    colnames(aux.IPI.c1)[i],".PNG", sep="")
  ggsave(p, file=file_name, width = 14, height = 10, units = "cm")
}

#Compute IPI Median
m.IPI.c1=lapply(aux.IPI.c1, median)
m.IPI.c2=lapply(aux.IPI.c2, median)
m.IPI.c3=lapply(aux.IPI.c3, median)
m.IPI.fastICA=lapply(aux.IPI.fastICA, median)
#Compute IPI Quantiles
q.IPI.c1=lapply(aux.IPI.c1, quantile,  probs=c(0.025,0.975))
q.IPI.c2=lapply(aux.IPI.c2, quantile,  probs=c(0.025,0.975))
q.IPI.c3=lapply(aux.IPI.c3, quantile,  probs=c(0.025,0.975))
q.IPI.fastICA=lapply(aux.IPI.fastICA, quantile,  probs=c(0.025,0.975))

###Bootstrap CI###
c1Results=data.frame()
c2Results=data.frame()
c3Results=data.frame()
fastICAResults=data.frame()
for (i in 1:3){
  #LR
  x1=aux.IPI.c1[,i]
  bootCI.IPI.c1=boot.ci(boot(x1, function(x,j) median(x[j]), R=10000))
  c1Results=rbind.data.frame(c1Results,data.frame(bootCI.IPI.c1$basic[4], median(x1), bootCI.IPI.c1$basic[5]))
  #Zhang
  x2=aux.IPI.c2[,i]
  bootCI.IPI.c2=boot.ci(boot(x2, function(x,j) median(x[j]), R=10000))
  c2Results=rbind.data.frame(c2Results,data.frame(bootCI.IPI.c2$basic[4], median(x2), bootCI.IPI.c2$basic[5]))
  #Lin
  x3=aux.IPI.c3[,i]
  bootCI.IPI.c3=boot.ci(boot(x3, function(x,j) median(x[j]), R=10000))
  c3Results=rbind.data.frame(c3Results,data.frame(bootCI.IPI.c3$basic[4], median(x3), bootCI.IPI.c3$basic[5]))
  #fastICA
  x4=aux.IPI.fastICA[,i]
  bootCI.IPI.fastICA=boot.ci(boot(x4, function(x,j) median(x[j]), R=10000))
  fastICAResults=rbind.data.frame(fastICAResults,data.frame(bootCI.IPI.fastICA$basic[4], median(x4), bootCI.IPI.fastICA$basic[5]))
}
colnames(c1Results)=c("c1.lowerCI", "c1.Median.IPI", "c1.upperCI")
colnames(c2Results)=c("c2.lowerCI", "c2.Median.IPI", "c2.upperCI")
colnames(c3Results)=c("c3.lowerCI", "c3.Median.IPI", "c3.upperCI")
colnames(fastICAResults)=c("fastICA.lowerCI", "fastICA.Median.IPI", "fastICA.upperCI")
IPI.Results=data.frame(c1Results, c2Results, c3Results, fastICAResults)
write.csv(IPI.Results,"C:/Users/sorel/Desktop/Paper ICA/Scripts/Simulation Results/IPI_Results_Sim4.csv", row.names = FALSE)













#TO DO
#1. Add Noice to the reference signals
#2. In case i dont know anything use the three classical load profiles curves. (generate functions with the approximate shape of the load curves)
#3. Pure signals of a particular type (data-driven approach). T°
#4. Data-Driven vs Domain Knowledge for the load profiles extraction (review domain driven load profile construction).
#5. review kernel-PCA (intractable with data growth)
#6. time aspect (how frequently are you gonna run these models)
#7. Common for walltime, number of flops IDEA: See how the algorithms scale 
#8. See what if's scenarios.
#9. Simulate the hierarchical data sets (how to specify the aggregation assume indepence and simulate)
#10. Dream-up arbitrarily scaleable situation (M references, N hidden signals)
#11. Construct the R package
#12. generate quantiles curves to see the CI for the simulations (Use commom random numbers)


