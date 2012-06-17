%Adaptive Monte Carlo
clear all, close all
tic

A=0.56; %constant in BE inequality
alpha=0.05; %uncertainty
alpha1=1-sqrt(1-alpha); %1-a1=sqrt(1-a)
zal2=norminv(1-alpha/2);
zal12=norminv(1-alpha1/2);
epsilon=0.01; %absolute error tolerance
n0=100;
L=1.2;
kurtmax=(n0-3)/(n0-1) + ((alpha1*n0)/(1-alpha1))*(1-1/L^2)^2;

disp(['The initial sample size for estimating variance = ' int2str(n0)])
disp(['The sample variance inflation factor = ' num2str(L^2)])
disp(['The maximum possible kurtosis = ' num2str(kurtmax)])



mu=1; %assumed mean
sig=1; %assumed standard deviation
c=0.01;

testfuntype='step';
if strcmp(testfuntype,'step')
    % Step test function 
    testfun = @(x,mu,sig,p) mu ...
        + (sig*sqrt((1-p)/p)).*(x<=p) ...
        - (sig*sqrt(p/(1-p))).*(x>p);
    kurt=1/(c*(1-c))-3;
    rhotype='unif';
elseif strcmp(testfuntype,'exp')
    % Exponential test function
    testfun = @(x,mu,sig,c) mu ...
        + (sig/sqrt(exp(c^2)-1))*(exp(c*x-c^2/2) -1);
    kurt=-3+3*exp(2*c^2)+2*exp(3*c^2)+exp(4*c^2);
    rhotype='norm';
end
disp(['The kurtosis of this example = ' num2str(kurt)])

wfN=zal2;
wfN1=zal12;
NN=ceil((wfN*sig/epsilon)^2);
disp(['Assuming normality and known variance, the correct sample size = ' int2str(NN)])

nrep=1000;
muhatBE=zeros(nrep,1); 
sig0=muhatBE; sig0up=muhatBE;
NNup=muhatBE; NB=muhatBE;
stddev=muhatBE;
for i=1:nrep  
    if strcmp(rhotype,'unif')
        x0=rand(n0,1);
    elseif strcmp(rhotype,'norm')
        x0=randn(n0,1);
    end
    y0=testfun(x0,mu,sig,c);
    sig0(i)=std(y0);
    sig0up(i)=sig0(i)*L;
    NNup(i)=ceil((wfN*sig0up(i)/epsilon)^2);
    if sig0up(i)==0;
        NB(i)=1;
    else
        NBEfun=@(logsqrtn) ...
            normcdf(-exp(logsqrtn).*epsilon/sig0up(i)) + ...
            A*kurtmax^(3/4).*exp(-logsqrtn).*(1+exp(logsqrtn).*epsilon/sig0up(i)).^(-3) - alpha1/2;
        logsqrtn=fzero(NBEfun,log(NNup(i)));
    NB(i)=ceil(exp(2*logsqrtn));
    end
    NB(i)=max(NB(i),n0);
    if strcmp(rhotype,'unif')
        x=rand(NB(i),1);
    elseif strcmp(rhotype,'norm')
        x=randn(NB(i),1);
    end
    y=testfun(x,mu,sig,c);
    muhatBE(i)=mean(y);
    stddev(i)=std(y);
end
pass=(stddev<=sig0up);
intol=abs(mu-muhatBE)<=epsilon;
probintol=mean(intol)
probpassintol=sum(intol&pass)/sum(pass)

toc
