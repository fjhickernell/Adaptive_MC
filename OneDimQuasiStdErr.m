% 1D Quasi-Standard Error experiments

%% Initialize
clear all, close all
Mavg=7;
mmax=10;
nvec=2.^(0:mmax)';
nmax=nvec(end)*Mavg;
fudge=1;
nplot=1000;

funname='sparsespect';
switch funname
    case 'hump'
        pow=3;
        testfun=@(x) (x.*(1-x)).^pow;
        syms xsym; fsym=testfun(xsym);
        exactInteg=double(int(fsym,0,1));
    case 'sparsespect'
        mMaxFreq=20;
        mskip=4;
        logFreqVec=(0:mskip:mMaxFreq)';
        freqVec=(2.^logFreqVec)*Mavg;
        testfun=@(x) 1+ ...
            sum(cos(((2*pi)*freqVec)*x).*repmat(1./(1+logFreqVec).^2,1,length(x)),1);
        exactInteg=1;
end

shift=rand(1,1)/nmax;
xnode=(0:nmax-1)./nmax+shift;
fnode=testfun(xnode);

%% Compute integral
appxInteg=zeros(mmax,1);
qse=appxInteg;
for m=0:mmax
    n=nvec(m+1);
    nskip=2^(mmax-m);
    fuse=fnode(1:nskip:nmax);
    fpart=reshape(fuse,Mavg,n);
    appxIntegPart=mean(fpart,2);
    appxInteg(m+1)=mean(appxIntegPart);
    qse(m+1)=std(appxIntegPart)/sqrt(Mavg);
end
trueErr=abs(exactInteg-appxInteg);
errEst=fudge*qse;

%% Plot
figure;
xplot=(0:nplot)/nplot;
fplot=testfun(xplot);
plot(xplot,fplot,'b-',...
    'linewidth',2)

figure;
loglog(nvec*Mavg,trueErr,'b.',...
    nvec*Mavg,errEst,'r.',...
    'markersize',20)



