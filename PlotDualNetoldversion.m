%Plot dual net
close all, clear all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

%% Set up wavenumbers
kexp=10;
kmax=2^kexp;
kvec=(0:kmax-1);
kvecplot=[1/2 (1:kmax-1)];
[kxx,kyy]=meshgrid(kvec);
[kxplot,kyplot]=meshgrid(kvecplot);

%% Digits of wavenumbers
kdig=zeros(kmax,kexp);
for m=1:kexp
    mm1=m-1;
    twomm1=2^mm1;
    kdig(twomm1+(1:twomm1),1:m)=...
        [kdig(1:twomm1,1:mm1) ones(twomm1,1)];
end

%% Plot wavenumbers
% figure;
% loglog(kxplot(:),kyplot(:),'b.')
% figure;
% plot(kxx(:),kyy(:),'b.')

%% Get net
mmax=10;
mvec=1:mmax;
kxdual=repmat(kxx,mmax);
kydual=repmat(kyy,mmax);
pSob=sobolset(2);
xSob=net(pSob,2^mmax);
xSobbase=xSob(2.^(0:mmax-1)+1,:);
xdig=bitget(repmat(floor(xSobbase(:,1)*kmax),1,kexp),repmat(1:kexp,mmax,1));
xdig=xdig(:,kexp:-1:1);
ydig=bitget(repmat(floor(xSobbase(:,2)*kmax),1,kexp),repmat(1:kexp,mmax,1));
ydig=ydig(:,kexp:-1:1);

%% Compute dot products
kdotxbase=mod(kdig*xdig',2);
kdotybase=mod(kdig*ydig',2);
kdotxvec=zeros(kmax*kmax,mmax);
for m=mvec
    kdotxtemp=repmat(kdotxbase(:,m),1,kmax);
    kdotytemp=repmat(kdotybase(:,m),1,kmax);
    kdotxmat=mod(kdotxtemp'+kdotytemp,2);
    kdotxvec(:,m)=kdotxmat(:);
    %keyboard
end
isdual=logical(cumprod(1-kdotxvec,2));

%% Plot dual lattices
r=3;
cut=[inf 1e-1 1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5 0];
ncut=length(cut)-1;
dotsize=[40 20 10 5 2 1 0.5 0.2 0.1];
for m=r+1:mmax
    figure;
    set(gca,'xtick',[0.5 1 10 100 1000],'ytick',[0.5 1 10 100 1000],...
        'xticklabel',[0 1 10 100 1000],'yticklabel',[0 1 10 100 1000],...
        'Xscale','Log','Yscale','Log')
    axis([0.5 kmax 0.5 kmax])
    hold on
    fulldual=isdual(:,m);
    fulldualp=fulldual; fulldualp(1)=false;
    partdual=isdual(:,m-r);
    diffdual=partdual&(~fulldual);
    kxdiff=kxplot(diffdual);
    kydiff=kyplot(diffdual);
    kxfull=kxplot(fulldualp);
    kyfull=kyplot(fulldualp);
    kdiffsize=1./(kxdiff.*kydiff);
    kfullsize=1./(kxfull.*kyfull);
    for j=1:ncut
        whdiff=(kdiffsize<cut(j) & kdiffsize>=cut(j+1));
        whfull=(kfullsize<cut(j) & kfullsize>=cut(j+1));          
        plot(kxdiff(whdiff),kydiff(whdiff),'b.',...,
            kxfull(whfull),kyfull(whfull),'r.',...       
            'markersize',dotsize(j))
        %keyboard
    end
    set(gca,'xtick',[0.5 1 10 100 1000],'ytick',[0.5 1 10 100 1000],...
        'xticklabel',[0 1 10 100 1000],'yticklabel',[0 1 10 100 1000])
    axis([0.5 kmax 0.5 kmax])
    xlabel('{\it k}_1')
    ylabel('{\it k}_2')
    eval(['print -depsc dualnet-' int2str(m) '.eps'])
    %keyboard
end

