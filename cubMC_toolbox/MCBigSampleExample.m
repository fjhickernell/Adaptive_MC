% Example to test MC for large numbers of samples
% Compute the price of a call option by Monte Carlo
% Rather than just using one long array, 
%   we try for loops to test the time required

clear all, format compact %garbage collection

r=0.01; %interest rate
sigma=0.5; % volatility
mu=r-sigma^2/2;
d=1; %number of time steps
S0=100; %initial stock price
K=100; %strike price
T=1; %time to expiry
dt=T/d; %delta t
sigdt=sigma*sqrt(dt);
trueprice=S0*normcdf((log(S0/K)+(mu+sigma^2)*T)/(sigma*sqrt(T))) ...
    -K*exp(-r*T)*normcdf((log(S0/K)+mu*T)/(sigma*sqrt(T)));

nsamp=1e8; %sample size

mvec=4:8; %log10 of sizes of arrays to try
nm=length(mvec);
time=zeros(nm,1);

disp(['True price = ' num2str(trueprice,'%3.2f')])
for jj=1:nm
    tic
    nsize=10^mvec(jj); %size of array
    narray=nsamp/nsize; %number of arrays
    meanpart=zeros(narray,1);
    sspart=zeros(narray,1);
    for ii=1:narray;
        x=randn(nsize,d); %call random numbers
        smat=S0*cumprod([ones(nsize,1) exp(mu*dt+sigdt*x)],2);
        payoff=max(smat(:,d+1)-K,0)*exp(-r*T);
        meanpart(ii)=mean(payoff);
        sspart(ii)=sum((payoff-meanpart(ii)).^2,1);
    end
    overallmean=mean(meanpart); %mean 
    ssmean=sum((meanpart-overallmean).^2,1);
    overallss=sum(sspart,1)+ssmean;
    overallstd=sqrt(overallss/(nsamp-1)); %standard deviation
    err=overallstd*(1.96/sqrt(nsamp)); %CLT error estimate
    time(jj)=toc;
    disp(['Price = ' num2str(overallmean,'%3.2f') ...
        ' +/- ' num2str(err,'%3.5f') ' in ' ...
        num2str(time(jj),'%3.3f') ' sec w/ ' ...
        int2str(narray) ' arrays of size ' int2str(nsize) ' x ' ...
        int2str(d)])
end

