%Adaptive MC comparison
%Some plots
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
clear all, close all

%% First a plot of sample sizes
amin=1e-5;
amax=0.5;
nal=100;
alphavec=exp(log(amin)+log(amax/amin)*(0:1:nal-1)/(nal-1));
za22=(norminv(alphavec/2)).^2;
rho=50;
tolovsig=0.01;
A=0.56;
n0=1000; %initial sample size
fudge=1.5; %fudge factor
kurtmaxfun=@(n,fudge,alpha) (n-3)/(n-1) + ((alpha*n)/(1-alpha))*(1 - 1/(fudge^2))^2;

nc=ceil(1./((tolovsig.^2).*alphavec));
nn=ceil(za22./(tolovsig.^2));
nb=zeros(1,nal);
nbnon=zeros(1,nal);
for i=1:nal
    alpha=alphavec(i);
    kurtmax=kurtmaxfun(n0,fudge,alpha);
    rho=kurtmax.^(3/4);
    NBEfun=@(logsqrtn) normcdf(-exp(logsqrtn).*tolovsig) + A*rho.*exp(-logsqrtn) - alpha/2;
    NBEnonfun=@(logsqrtn) normcdf(-exp(logsqrtn).*tolovsig) + ...
        A*rho.*exp(-logsqrtn)./(1+exp(logsqrtn).*tolovsig).^3 - alpha/2;
    logsqrtn=fzero(NBEfun,log(sqrt(nn(i))));
    nb(i)=ceil(exp(2*logsqrtn));
    logsqrtn=fzero(NBEnonfun,log(sqrt(nn(i))));
    nbnon(i)=ceil(exp(2*logsqrtn));
end
amaxexp=ceil(log10(amax)); aminexp=floor(log10(amin));
Nmaxexp=ceil(log10(max([nc nn nbnon]))); Nminexp=floor(log10(min([nc nn nbnon])));

figure;
h=loglog(alphavec,nc,'b-',alphavec,nn,'g-',alphavec,nbnon,'r-','linewidth',2);
xlabel('{\it \alpha}')
set(gca,'Xtick',10.^(aminexp:amaxexp),'Ytick',10.^(Nminexp:Nmaxexp))
axis(10.^[aminexp amaxexp Nminexp Nmaxexp])
legend(h,{'{\it N_C}','{\it N_G}','{\it N_B}'})
print -depsc alphacompare.eps


figure;
h=loglog(alphavec,nn,'g-',alphavec,nbnon,'r-','linewidth',2);
xlabel('{\it \alpha}')
set(gca,'Xtick',10.^(aminexp:amaxexp),'Ytick',10.^(Nminexp:Nmaxexp))
Nmaxexp=ceil(log10(max([nn nbnon]))); Nminexp=floor(log10(min([nn nbnon])));
axis(10.^[aminexp amaxexp Nminexp Nmaxexp])
legend(h,{'{\it N_G}','{\it N_{B}}'})
print -depsc alphacompareb.eps

L=2;
alpha=0.05;
M=(L/(L-1))*((1-alpha)/alpha)
nmin=(1+sqrt(1+8*M))/2

%% Then a plot of kutosis max
nvec=[30 100 1000 1e4 1e5]';
nn=length(nvec);
L=1.5;
amin=1e-3;
amax=0.5;
nal=100;
alphavec=exp(log(amin)+log(amax/amin)*(0:1:nal-1)/(nal-1));
temp1=repmat((nvec-3)./(nvec-1),1,nal);
temp2=nvec*(alphavec*(1-1/L^2)^2./(1-alphavec));
kmax=temp1+temp2;

figure;
h=loglog(alphavec,kmax,'-','linewidth',2);
legend(h,{['{\it n}_\sigma = ' num2str(nvec(1))],...
    ['{\it n}_\sigma = ' num2str(nvec(2))],...
    ['{\it n}_\sigma = ' num2str(nvec(3))],...
    ['{\it n}_\sigma = ' num2str(nvec(4))],...
    ['{\it n}_\sigma = ' num2str(nvec(5))]},...
    'location','northwest');
%xlabel('{\it \alpha}')
text('Position',[2e-2,3e-1],'Interpreter','Latex',...
    'String','$$\tilde{\alpha}$$')
ylabel('{\it \kappa}_{max}')
axis(10.^[-3 0 0 5])
print -depsc kurtmaxfig.eps
alphagood=0.1;
temp3=nvec*(alphagood*(1-1/L^2)^2./(1-alphagood));
kmaxgood=temp3+(nvec-3)./(nvec-1);
disp(['For alpha = ' num2str(alphagood)])
disp(['n = ' int2str(nvec')])
disp(['kmax = ' num2str(kmaxgood')])

%% A plot on the upper bound of the probability that v_n >= hsigma^2
n0vec=[30 100 1000 1e4]';
nvec=2.^(10:20);
nn0=length(n0vec);
nn=length(nvec);
n0mat=repmat(n0vec,1,nn);
nmat=repmat(nvec,nn0,1);
L=1.5;
alpha=0.05;
al1=1 - sqrt(1-alpha);
temp1=1+(n0mat./(L^4*nmat));
temp2=al1/(1-al1);
temp3=2*(nmat-n0mat)./(nmat.*(n0mat-1).*(nmat-1).*(L^2-1)^2);
temp4=temp1*temp2 - temp3;
prob=1./(1+1./temp4);

figure;
h=semilogx(nvec,prob,'-','linewidth',2);
legend(h,{['{\it n}_\sigma = ' num2str(n0vec(1))],...
    ['{\it n}_\sigma = ' num2str(n0vec(2))],...
    ['{\it n}_\sigma = ' num2str(n0vec(3))],...
    ['{\it n}_\sigma = ' num2str(n0vec(4))]},'location','northeast')
xlabel('{\it n}')
ylabel('Probability')
axis([10^3 10^6 0 0.1])
print -depsc probfailfig.eps


%% Then a plot of the amplification factor
% nvec=[30 100 1000 1e4]';
% kurt=5;
% fac1=(kurt-(nvec-3)./(nvec-1))./nvec;
% alphaasym=1./(1+1./fac1);
% size(alphavec)
% alphavec=sort([alphavec alphaasym'*1.001]);
% size(alphavec)
% ampfac=1./(1 - sqrt(fac1*((1-alphavec)./alphavec)));
% whpos=ampfac>0;
% 
% figure;
% h=loglog(alphavec(whpos(1,:)),ampfac(1,whpos(1,:)),'b-',...
%     alphavec(whpos(2,:)),ampfac(2,whpos(2,:)),'g-',...
%     alphavec(whpos(3,:)),ampfac(3,whpos(3,:)),'r-',...
%     alphavec(whpos(4,:)),ampfac(4,whpos(4,:)),'k-','linewidth',2);
% legend(h,{['{\it n} = ' num2str(nvec(1))],...
%     ['{\it n} = ' num2str(nvec(2))],...
%     ['{\it n} = ' num2str(nvec(3))],...
%     ['{\it n} = ' num2str(nvec(4))]},'location','northwest')
% xlabel('{\it \alpha}')
% ylabel('Amplification Factor')
% axis(10.^[-4 0 0 3])
% print -depsc varianceampfac.eps
