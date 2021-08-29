library(MASS) 
library(fda)
library(refund)

fun_gendata=function(n,pn,ppt,sig,ex,mu)
{
  # n is the number of sample size
  # k is the number of FPC scores
  # pn is the number of scalars
  # ppt is the dimension of functional data
  K=50
  pn4=pn+K;   
  sigma=matrix(0,pn4,pn4);  
  for (i in 1: K)
  {sigma[i,i]=16*i^(-1)}
  
  for(i in 1:K){
    for(j in (K+1):pn4){
      sigma[i,j]=0.2^(abs(i-j+K)+1);  
      sigma[j,i]=sigma[i,j];}}
  
  for(i in (K+1):pn4){
    for(j in (K+1):pn4){
      sigma[i,j]=0.5^(abs(i-j));}  # Z: AR(0.5) 
  }
  
  txi_z=mvrnorm(n=n,rep(0,pn4),sigma) 
  xi=txi_z[,1:K]; 
  z=txi_z[,(K+1):pn4]  
  ###########################################################################

   phi=function(s,k){
      if(k%%2==0)  return(sqrt(2)*sin((k-1)*pi*s)/sqrt(ppt))  else  return(sqrt(2)*cos(k*pi*s)/sqrt(ppt));}
  
  kk=matrix(rep(c(1:K),ppt),ppt,K,byrow=TRUE);  
  tt=matrix(rep(seq(0,1,length=ppt),K),ppt,K,byrow=FALSE); 
  phi_kt=matrix(mapply(function(s,k)phi(s,k),tt,kk),ppt,K);  # ppt*k   
  x_true=xi%*%t(phi_kt)  #n*ppt
  w=x_true+matrix(rnorm(n*ppt,0,0.5),n,ppt,byrow = TRUE); # n*ppt
  #############################################################################
  b=rep(0,K)
for (j in 1:K)
{
b[j]=4*(-1)^(j+1)*j^(-2)
}
b[1]=0.3
 # b=c(2,1,-1,0.5) # 1*k
  beta_true=t(matrix(b))%*%t(phi_kt)  #1*ppt
  gamma_true=rep(1,pn)  # 1*pn vector which is used for scalar variates
  epsilon=rnorm(n,0,sig)
  y=xi%*%matrix(b)+z%*%matrix(gamma_true)+matrix(epsilon)
  ###############################################################################
  if (ex==1){
             u=rbinom(n, 1,0.5)
             u=2*u-1
             mu_true=u%*%t(mu)
     ys=y+mu_true
    }
   if (ex==2){
             u=sample(c(-1:1),n,replace = TRUE, prob = c(1/3,1/3,1/3))
             mu_true=u%*%t(mu)
     ys=y+mu_true
    }
   if (ex==3){

             mu_true=rep(1,n)%*%t(mu)
     ys=y+mu_true
    }
  list(y=ys,z=z,w=w,gamma_true=gamma_true,beta_true=beta_true,mu_true=mu_true,u=u)
}
#################################################################

initial<-function(z,x,y,n,lamin,al)
{

dx=ncol(x)
dz=ncol(z)

 Qz=diag(n)-z%*%solve(t(z)%*%z)%*%t(z)
 
DEL=matrix(0,nrow=n,ncol=(n*(n-1)/2))  ##for e_i-e_j
e=diag(n)
XD=t(x[1,])
for (i in 1:(n-1)){
ind=(i-1)*n-i*(i-1)/2
DEL[,(ind+1):(ind+n-i)]=e[,i]-e[,(i+1):n]   ##for n beta(p*1)
XD=bdiag(XD,t(x[(i+1),]))
}
DEL=t(DEL)
AA=t(DEL)%*%DEL
AA=kronecker(AA,diag(dx))

pr=t(XD)%*%Qz%*%XD+lamin*AA
pr=solve(pr)
betain=pr%*%t(XD)%*%Qz%*%y
etain=solve(t(z)%*%z)%*%t(z)%*%(y-XD%*%betain)

betain=matrix(betain,nrow=dx,ncol=n)
betam=apply(betain,2,median)
#betam=betain[2,]
betaind=cbind(betam,seq(1,n))
betaind=betaind[order(betaind[,1]),]
Kstar=floor(0.5*n^(1/2))
GA=round(n/Kstar)

W=matrix(0,nrow=n,ncol=Kstar)

for(i in 1:(Kstar-1)){
iid=betaind[((i-1)*GA+1):(i*GA),2]
W[iid,i]=1
}
iid=betaind[((Kstar-1)*GA+1):n,2]
W[iid,Kstar]=1

W=kronecker(W,diag(dx))
XW=XD%*%W
XW=as.matrix(XW)
ZWX=cbind(z,XW)
#thetain=solve(t(ZWX)%*%ZWX+0.01*diag(ncol(ZWX)))%*%t(ZWX)%*%y
thetain=ginv(t(ZWX)%*%ZWX)%*%t(ZWX)%*%y
etain=thetain[1:dz]
alphain=thetain[(dz+1):length(thetain)]
betain=W%*%alphain

return(list(Qz,XD,AA,DEL,betain,etain))
  }

#########################################################################
lamoptimal<-function(z,x,y,n,lam,varth,gam,al,betaini)
{
      betain=betaini
	lamU=4
	lamL=1
	#lamU=2
	#lamL=0.2
	lams=seq(lamL,lamU,length=10)
	ss=length(lams)
	BIC=seq(-1,-1,length=ss)
	deg=seq(-1,-1,length=ss)
      	Qnn=seq(-1,-1,length=ss)
     betas=matrix(0,nrow=n*ncol(x),ncol=ss)
	for(i in 1:ss){
		lam=lams[i]
		result=estimation(z,x,y,n,lam,varth,gam,al,betain,etain,Qz,XD,AA,DEL)
            etahat=result[[3]]
            betahat=result[[2]]
            betain=result[[2]]
            betain=as.vector(betain)
            betas[,i]=betain
            alphahat=result[[5]]
            W=result[[6]]
            sighat=result[[7]]
            Khat=result[[1]]
            XD=result[[8]]
            df=Khat*ncol(x)+ncol(z)
            Qn=t(y-z%*%etahat-XD%*%W%*%alphahat)%*%(y-z%*%etahat-XD%*%W%*%alphahat)/n
            Qn=Qn[1,1]
            Qnn[i]=Qn
		BIC[i]=log(Qn)+10*log(log(n*ncol(x)+ncol(z)))*log(n)*df/n

	}
	

      BIC=BIC[-1]
	lam=lams[(max(which(BIC==min(BIC)))+1)]
	ind=max(which(BIC==min(BIC)))+1
    betain=betas[,ind]
     
	
	    return(list(lam,BIC,betain))
  }
###########################################################################################



estimation<-function(z,x,y,n,lam,varth,gam,al,betain,etain,Qz,XD,AA,DEL)
{
 dx=ncol(x)
 dz=ncol(z)
 betaold=betain
 etaold=etain
 upsilonold=matrix(0,nrow=(n*(n-1)/2),ncol=dx)
 upsilonold=t(upsilonold)
 #upsilonold=as.vector(upsilonold)
 betat=matrix(betaold,nrow=dx,ncol=n)
 betat=t(betat)
 deltaold=DEL%*%betat 
 deltaold=t(deltaold)
 #deltaold=as.vector(deltaold)

  B=t(XD)%*%Qz%*%XD+varth*AA
  B=solve(B)
  ep=10^(-3)
  step=0
  lam1=lam/varth

  rr=10

  while(rr>ep)
  {
  step=step+1 
  du=deltaold-upsilonold*varth^(-1)
  dD=du%*%DEL
  vdD=as.vector(dD)
  betanew=t(XD)%*%Qz%*%y+varth*vdD
  betanew=B%*%betanew
  etanew=solve(t(z)%*%z)%*%t(z)%*%(y-XD%*%betanew)
  
 betat=matrix(betanew,nrow=dx,ncol=n)
 betat=t(betat)
 deltat=DEL%*%betat 
 deltat=t(deltat)
 zeta=deltat+varth^(-1)*upsilonold
 zeta=t(zeta)

 zetanorm=apply(zeta^2,1,sum)
  zetanorm=sqrt(zetanorm)
  SS=1-lam1/zetanorm
  S=(((SS>0)*SS)%*%t(seq(1,1,length=dx)))*zeta
  thr1=(zetanorm>(gam*lam))%*%t(seq(1,1,length=dx))
  thr2=(zetanorm<=(gam*lam))%*%t(seq(1,1,length=dx))

  lam2=lam1*gam/(gam-1)
  SS2=1-lam2/zetanorm
  S2=(((SS2>0)*SS2)%*%t(seq(1,1,length=dx)))*zeta
  thr3=(zetanorm<=(lam+lam1))%*%t(seq(1,1,length=dx))
  thr4=((zetanorm>(lam+lam1))*(zetanorm<=(gam*lam)))%*%t(seq(1,1,length=dx))

  if(al==1){
  deltanew=S
  }else if(al==2){
  deltanew=zeta*thr1+(S/(1-(varth*gam)^(-1)))*thr2
  }else{
  deltanew=zeta*thr1+S*thr3+(S2/(1-((gam-1)*varth)^(-1)))*thr4
  }
  deltanew=t(deltanew)
  #deltanew=as.vector(deltanew)
  upsilonnew=upsilonold+varth*(deltat-deltanew)
  r=deltat-deltanew
  r=as.vector(r)
  rr=sqrt(t(r)%*%r)

  betaold=betanew
  etaold=etanew
  deltaold=deltanew
  upsilonold=upsilonnew

}


bb=matrix(betanew,nrow=dx,ncol=n)
bb=t(bb)

deltam=t(deltanew)
deltan=apply(abs(deltam),1,sum)

seq=1:n
group=matrix(0,nrow=n,ncol=n)
alphaold=matrix(0,nrow=n,ncol=dx)
Wtilda=matrix(0,nrow=n,ncol=n)
K=1

while(length(seq)>0)
{
i=seq[1]
ind=(i-1)*n-i*(i-1)/2
id=which(deltan[(ind+1):(ind+n-i)]==0)
id=c(i,i+id)
group[1:length(id),K]=id
if(length(id)>1){
bbm=bb[id,]
bbm=as.matrix(bbm)
alphaold[K,]=apply(bbm,2,mean)
}else{
alphaold[K,]=bb[id,]
}
Wtilda[id,K]=1
seq=seq [! seq %in% id]
K=K+1
}

K=K-1


alphaold=alphaold[1:K,]
Wtilda=Wtilda[,1:K]
Wtilda=as.matrix(Wtilda)
W=kronecker(Wtilda,diag(dx))
alphaold=t(alphaold)
alphaold=as.vector(alphaold)

sigg=sqrt(t(y-z%*%etaold-XD%*%W%*%alphaold)%*%(y-z%*%etaold-XD%*%W%*%alphaold)/(n-dz-K*dx))
a=a0=matrix(rep(1,n*n),n,n)
Phat=P0=matrix(0,n,n)
for (i in 1:(n-1))
{for (j in (i+1):n)
{a[i,j]=sum((W[i,]-W[j,])^2)
a0[i,j]=mu_true[i]-mu_true[j]}
Phat[i,(which(a[i,]==0))]=1
P0[i,(which(a0[i,]==0))]=1
}
RI=(length(which((P0-Phat)==0))-n-n*(n-1)/2)/(n*(n-1)/2)

return(list(K,betaold,etaold,group,alphaold,W,sigg,XD,Wtilda,RI))
  }
  
###########################################################



n=400
varth=1
gam=3
pn=5
ppt=100
sig=1
nsim=500#200

ex=1
mu=3
al=3 #2 MCP 1 LASSO 3 SCAD
ph=1
MSEgamma=MSEmu=MSEbeta=MSEgamma_lse=MSEmu_lse=MSEbeta_lse=MSEgamma_opt=MSEmu_opt=MSEbeta_opt=Ksim=RI=lamb=rep(0,nsim)
###################################################

for (sim in 1:nsim)
{
######################  generate data
DAT=fun_gendata(n,pn,ppt,sig,ex,mu)
y=DAT$y
z=DAT$z
w=DAT$w
u=DAT$u
gamma_true=DAT$gamma_true
beta_true=DAT$beta_true
mu_true=DAT$mu_true
##########################

###########################
####  fpca dimension reduction
pcs=fpca.face(Y = w, argvals = seq(0,1,len=ppt),pve = 0.95)
sn=pcs$npc  #number of pc
xi=pcs$scores  #n*sn
phihat=pcs$efunctions  #ppt*sn
######new z in huang
znew=cbind(xi,z)
x=matrix(rep(1,n),n,1)
#################


###### olse
estimate_lse=lm(y~znew)
mu_lse=estimate_lse$coefficients[1]
w_lse=estimate_lse$coefficients[c(2:(sn+1))]
z_lse=estimate_lse$coefficients[c((sn+2):length(estimate_lse$coefficients))]

betahatf_lse=t(w_lse)%*%t(phihat)
muhat_lse=rep(mu_lse,n)
##########################################################
#####  orecle estimate
u_true=unique(u)
len_u=length(u_true)
ww=matrix(0,n,len_u)
for (uu in 1:len_u)
{
ww[which(u==u_true[uu]),uu]=1
}
znew_opt=cbind(ww,znew)

estimate_opt=lm(y~znew_opt-1)
mu_opt=estimate_opt$coefficients[c(1:(len_u))]
w_opt=estimate_opt$coefficients[c((len_u+1):(len_u+sn))]
z_opt=estimate_opt$coefficients[c((len_u+1+sn):length(estimate_opt$coefficients))]
betahatf_opt=t(w_opt)%*%t(phihat)

muhat_opt=matrix(0,n,1)
for (ii in 1:n)
{
muhat_opt[ii,]=mu_opt[which(ww[ii,]==1)]
}


####################### subgroup analysis
lamin=0.001
results=initial(znew,x,y,n,lamin,al)
Qz=results[[1]]
XD=results[[2]]
AA=results[[3]]
DEL=results[[4]]
betaini=results[[5]]
etain=results[[6]]

lamopt=lamoptimal(znew,x,y,n,lam,varth,gam,al,betaini)
lam=lamopt[[1]]
betain=lamopt[[3]]

result=estimation(znew,x,y,n,lam,varth,gam,al,betain,etain,Qz,XD,AA,DEL)

etahat=result[[3]]
etahat=as.vector(etahat)
etahat_xi=etahat[c(1:sn)]
betahatf=t(etahat_xi)%*%t(phihat)

etahat_z=etahat[c((sn+1):length(etahat))]

betahat=result[[2]]
betahat=as.vector(betahat)
##########################################################

RI[sim]=result[[10]]
Ksim[sim]=result[[1]]
MSEbeta[sim]=sqrt((betahatf-beta_true)%*%t(betahatf-beta_true)/ppt)
MSEgamma[sim]=sqrt(t(etahat_z-gamma_true)%*%(etahat_z-gamma_true)/pn)
MSEmu[sim]=sqrt(t(betahat-mu_true)%*%(betahat-mu_true)/n)
lamb[sim]=lam

MSEbeta_lse[sim]=sqrt((betahatf_lse-beta_true)%*%t(betahatf_lse-beta_true)/ppt)
MSEgamma_lse[sim]=sqrt(t(z_lse-gamma_true)%*%(z_lse-gamma_true)/pn)
MSEmu_lse[sim]=sqrt(t(muhat_lse-mu_true)%*%(muhat_lse-mu_true)/n)

MSEbeta_opt[sim]=sqrt((betahatf_opt-beta_true)%*%t(betahatf_opt-beta_true)/ppt)
MSEgamma_opt[sim]=sqrt(t(z_opt-gamma_true)%*%(z_opt-gamma_true)/pn)
MSEmu_opt[sim]=sqrt(t(muhat_opt-mu_true)%*%(muhat_opt-mu_true)/n)
}
AAA=list(Ksim,RI,lamb,MSEmu,MSEbeta,MSEgamma,MSEmu_lse,MSEbeta_lse,MSEgamma_lse,MSEmu_opt,MSEbeta_opt,MSEgamma_opt)
meank=mean(Ksim)
sdk=sd(Ksim)
mediank=median(Ksim)

meanri=mean(RI)
sdri=sd(RI)

meanmsemu=mean(MSEmu)
sdmsemu=sd(MSEmu)
meanmsebeta=mean(MSEbeta)
sdmsebeta=sd(MSEbeta)
meanmsegamma=mean(MSEgamma)
sdmsegamma=sd(MSEgamma)

meanmsemu_lse=mean(MSEmu_lse)
sdmsemu_lse=sd(MSEmu)
meanmsebeta_lse=mean(MSEbeta_lse)
sdmsebeta_lse=sd(MSEbeta_lse)
meanmsegamma_lse=mean(MSEgamma_lse)
sdmsegamma_lse=sd(MSEgamma_lse)

meanmsemu_opt=mean(MSEmu_opt)
sdmsemu_opt=sd(MSEmu_opt)
meanmsebeta_opt=mean(MSEbeta_opt)
sdmsebeta_opt=sd(MSEbeta_opt)
meanmsegamma_opt=mean(MSEgamma_opt)
sdmsegamma_opt=sd(MSEgamma_opt)
BBB_mean=list(meank=meank,mediank=mediank,meanri=meanri,meanmsemu=meanmsemu,meanmsebeta=meanmsebeta,meanmsegamma=meanmsegamma,meanmsemu_lse=meanmsemu_lse,meanmsebeta_lse=meanmsebeta_lse,meanmsegamma_lse=meanmsegamma_lse,meanmsemu_opt=meanmsemu_opt,meanmsebeta_opt=meanmsebeta_opt,meanmsegamma_opt=meanmsegamma_opt)
BBB_sd=list(sdk=sdk,sdri=sdri,sdmsemu=sdmsemu,sdmsebeta=sdmsebeta,sdmsegamma=sdmsegamma,sdmsemu_lse=sdmsemu_lse,sdmsebeta_lse=sdmsebeta_lse,sdmsegamma_lse=sdmsegamma_lse,sdmsemu_opt=sdmsemu_opt,sdmsebeta_opt=sdmsebeta_opt,sdmsegamma_opt=sdmsegamma_opt)
BBB_mean
BBB_sd
setwd("C:/Users/86155/Desktop/subgroup_simu")
write.table(BBB_mean,"n400BBB_mean_ex1sig1mu3.csv",sep=",")
write.table(BBB_sd,"n400BBB_sd_ex1sig1mu3.csv",sep=",")
write.table(AAA,"n400AAA_ex1sig1mu3.csv",sep=",")

#############################################################