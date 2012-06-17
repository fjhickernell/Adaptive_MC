%Adaptive Monte Carlo
clear all, close all
tic

alpha=0.05; %uncertainty
alpha1=1-sqrt(1-alpha); %1-a1=sqrt(1-a)
zal2=norminv(1-alpha/2);
zal12=norminv(1-alpha1/2);
epsilon=0.01; %absolute error tolerance
mu=1; %assumed mean
sig=1; %assumed standard deviation
p=0.5;
n0=1000;
disp(['The initial sample size for estimating variance = ' int2str(n0)])
kurtposs=n0*alpha1/(1-alpha1) + 2*(n0-3)/(n0-1);
disp(['The maximum possible maximum kurtosis = ' num2str(kurtposs)])
kurtmax=min(0.9*kurtposs,2);
disp(['The maximum kurtosis assumed = ' num2str(kurtmax)])

testfun = @(x,mu,sig,p) mu ...
    + (sig*sqrt((1-p)/p)).*(x<=p) ...
    - (sig*sqrt(p/(1-p))).*(x>p);
kurt=1/(p*(1-p))-3;
disp(['The kurtosis of this example = ' num2str(kurt)])

wfN=zal2;
wfN1=zal12;
NN=ceil((wfN*sig/epsilon)^2);
disp(['Assuming normality and known variance, the correct sample size = ' int2str(NN)])

wfC=1/sqrt(alpha);
wfC1=1/sqrt(alpha1);
NC=ceil((wfC*sig/epsilon)^2);
disp(['Using Chebyshev' char(39) ' innequality and known variance, the correct sample size = ' int2str(NC)])

nrep=100;
muhatN=zeros(nrep,1); muhatNappx=muhatN; muhatNup=muhatN;
muhatC=muhatN; muhatCappx=muhatN; muhatCup=muhatN;
ciN=zeros(nrep,2); ciNappx=ciN; ciNup=ciN;
ciC=ciN; ciCappx=ciN; ciCup=ciN;
sig0=muhatN; sig0up=muhatN;
NNappx=muhatN; NNup=muhatN; NCappx=muhatN; NCup=muhatN;
upfac=1/sqrt(1-sqrt((kurtmax-(n0-3)/(n0-1))*(1-alpha1)/(alpha1*n0)));
disp(['The factor multiplying the sample variance to get the upper bound = ' num2str(upfac)])
for i=1:nrep  
    x0=rand(n0,1);
    y0=testfun(x0,mu,sig,p);
    sig0(i)=std(y0);
    sig0up(i)=sig0(i)*upfac;
    NNappx(i)=ceil((wfN*sig0(i)/epsilon)^2);
    NCappx(i)=ceil((wfC*sig0(i)/epsilon)^2);
    NNup(i)=ceil((wfN*sig0up(i)/epsilon)^2);
    NCup(i)=ceil((wfC*sig0up(i)/epsilon)^2);
    [muhatN(i),ciN(i,:)]=MCrun(testfun,NN,mu,sig,sig,p,wfN);
    [muhatNappx(i),ciNappx(i,:)]=MCrun(testfun,NNappx(i),mu,sig,sig0(i),p,wfN);
    [muhatNup(i),ciNup(i,:)]=MCrun(testfun,NNup(i),mu,sig,sig0(i),p,wfN1);
    [muhatC(i),ciC(i,:)]=MCrun(testfun,NC,mu,sig,sig,p,wfC);
    [muhatCappx(i),ciCappx(i,:)]=MCrun(testfun,NCappx(i),mu,sig,sig0(i),p,wfC);
    [muhatCup(i),ciCup(i,:)]=MCrun(testfun,NCup(i),mu,sig,sig0(i),p,wfC1);
    
end
probinciN=mean((ciN(:,1) <= mu) & (ciN(:,2) >= mu))
probinciNappx=mean((ciNappx(:,1) <= mu) & (ciNappx(:,2) >= mu))
probinciNup=mean((ciNup(:,1) <= mu) & (ciNup(:,2) >= mu))
probinciC=mean((ciC(:,1) <= mu) & (ciC(:,2) >= mu))
probinciCup=mean((ciCup(:,1) <= mu) & (ciCup(:,2) >= mu))

toc
