%Adaptive Monte Carlo
clear all, close all
tic

A=0.56;
alpha=0.05; %uncertainty
alpha1=1-sqrt(1-alpha); %1-a1=sqrt(1-a)
zal2=norminv(1-alpha/2);
zal12=norminv(1-alpha1/2);
epsilon=0.01; %absolute error tolerance
mu=1; %assumed mean
sig=1; %assumed standard deviation
n0=100;
kurtposs=n0*alpha1/(1-alpha1) + (n0-3)/(n0-1);
kurtmax=min(0.9*kurtposs,5);
nrep=1000;
pvec=[0.1 1 2];
np=length(pvec);
kurtvec=-3+exp(2*pvec.*pvec).*(3 + exp(pvec.*pvec).*(2 + exp(pvec.*pvec).*4));

testfun = @(x,mu,sig,p) mu - sig/(sqrt(exp(p.*p)-1)) + (sig/(sqrt(exp(p.*p)-1))).*exp(p*x-p.*p/2);

wfN=zal2;
wfN1=zal12;
NN=ceil((wfN*sig/epsilon)^2);
wfC=1/sqrt(alpha);
wfC1=1/sqrt(alpha1);
NC=ceil((wfC*sig/epsilon)^2);
disp(' ')

disp(['Uncertainty = ' num2str(alpha)])
disp(['Error tolerance required = ' num2str(epsilon)])
disp(['True mean = ' num2str(mu)])
disp(['True standard deviation = ' num2str(sig)])
disp(['Initial sample size for estimating variance = ' int2str(n0)])
disp(['Maximum possible kurtosis = ' num2str(kurtposs)])
disp(['Maximum kurtosis assumed = ' num2str(kurtmax)])
disp(['Assuming normality and known variance, the correct sample size = ' int2str(NN)])
disp(['Using Chebyshev' char(39) ' innequality and known variance, the correct sample size = ' int2str(NC)])


muhatBC=zeros(nrep,np);
sig0=zeros(nrep,1); sig0up=muhatBC;
errBC=muhatBC;
NNup=sig0; NCup=sig0; NBup=sig0; NCB=sig0;
probsucc=zeros(1,np);
upfac=1/sqrt(1-sqrt((kurtmax-(n0-3)/(n0-1))*(1-alpha1)/(alpha1*n0)));
disp(['The factor multiplying the sample variance to get the upper bound = ' num2str(upfac)])

for ii=1:np
    p=pvec(ii);
    kurt=kurtvec(ii);  
    disp(' ')
    disp(['True kurtosis = ' num2str(kurt)])
    for i=1:nrep  
        x0=rand(n0,1);
        y0=testfun(x0,mu,sig,p);
        sig0(i)=std(y0);
        sig0up(i,ii)=sig0(i)*upfac;
        NNup(i)=max(ceil((wfN1*sig0up(i,ii)/epsilon)^2),1);
        NCup(i)=max(ceil((wfC1*sig0up(i,ii)/epsilon)^2),1);
        if sig0up(i)==0;
            NBup(i)=1;
        else
            NBEfun=@(logsqrtn) normcdf(-exp(logsqrtn).*epsilon/sig0up(i,ii)) + A*kurtmax^(3/4).*exp(-logsqrtn) - alpha1/2;
            logsqrtn=fzero(NBEfun,log(NNup(i)));
            NBup(i)=ceil(exp(2*logsqrtn));
        end
        NCB(i)=min(NCup(i),NBup(i));
        x=rand(NCB(i),1);
        y=testfun(x,mu,sig,p);
        muhatBC(i,ii)=mean(y);
    end
errBC(:,ii)=abs(mu-muhatBC(:,ii));
probsucc(ii)=mean(errBC(:,ii) <= epsilon);
sortN=sort(NCB);
probcost=sortN(ceil((1-alpha)*nrep));
disp(['Probability of estimate within tolerance = ' num2str(probsucc(ii)*100) '%'])
disp(['With a ' num2str((1-alpha)*100) '% worst cost = ' int2str(probcost)])
end

save AdaptiveMCExpNormOut.mat

toc
