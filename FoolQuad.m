%Fool MATLAB's quad routine with a carefully chosen function

%% Initialize
clear all, close all %garbage collection
format compact, format long e %eliminate blank lines
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20) %make font labels large enough

nplot=500; %number of points to plot
tol=1e-14;

%% Set up function to fool quad
h=0.13579; %secret constant in quad

%A function that is constant at the ends and oscillatory in the middle
testfun=@(x) 1+cos(8*pi*min(max((x-2*h)/(1-4*h),0),1));

%Exact integral of this function
tic;
syms xsym hsym; hsym=h;
fsym=1+cos(8*pi*(xsym-2*hsym)/(1-4*hsym)); %expression in the middle
exactInteg=double(int(fsym,2*hsym,1-2*hsym)... %middle part
    +8*hsym); %plus end part
elapsedexacttime=toc;

%Use quad to find the answer
tic
[quadInteg,quadn]=quad(testfun,0,1,tol);
elapsedquadtime=toc;

%Use quadgk to find the answer
tic
[quadgkInteg,quadgkerr]=quadgk(testfun,0,1,'AbsTol',tol,'RelTol',0);
elapsedquadgktime=toc;

%Use quadl to find the answer
tic
[quadlInteg,quadln]=quadl(testfun,0,1,tol);
elapsedquadltime=toc;

%Use cubMC with iid sampling to find the answer
tic
paramiid.tol=1e-3;
paramiid.sample='iid';
[cubMCiidInteg,paramiid]=cubMC(testfun,[0 1],paramiid);
elapsedcubMCiidtime=paramiid.time;

%Use cubMC with Sobol sampling to find the answer
tic
paramSobol.tol=1e-3;
paramSobol.sample='Sobol';
[cubMCSobolInteg,paramSobol]=cubMC(testfun,[0 1],paramSobol);
elapsedcubMCSoboltime=paramSobol.time;

disp('For the function with well-placed wiggles')
disp(['Exact integral = ' num2str(exactInteg,'%10.15g')])
disp(['     which took ' num2str(elapsedexacttime) ' seconds to compute'])
disp(['Approximate integral using quad = ' num2str(quadInteg,'%10.15g') ...
    ' +/- ' num2str(tol,'%3.1g')])
disp(['     which took ' num2str(elapsedquadtime) ' seconds and'])
disp(['     ' num2str(quadn) ' function evaluations to compute'])
disp(['Approximate integral using quadgk = ' num2str(quadgkInteg,'%10.15g') ...
    ' +/- ' num2str(quadgkerr,'%3.1g')])
disp(['     which took ' num2str(elapsedquadgktime) ' seconds to compute'])
disp(['Approximate integral using quadl = ' num2str(quadlInteg,'%10.15g') ...
    ' +/- ' num2str(tol,'%3.1g')])
disp(['     which took ' num2str(elapsedquadltime) ' seconds and'])
disp(['     ' num2str(quadln) ' function evaluations to compute'])
disp(['Approximate integral using cubMC with iid sampling = '...
    num2str(cubMCiidInteg,'%10.15g') ...
    ' +/- ' num2str(paramiid.tol,'%3.1g')])
disp(['     which took ' num2str(elapsedcubMCiidtime) ' seconds and'])
disp(['     ' num2str(paramiid.n) ' function evaluations to compute'])
disp(['Approximate integral using cubMC with Sobol sampling = '...
    num2str(cubMCSobolInteg,'%10.15g') ...
    ' +/- ' num2str(paramSobol.tol,'%3.1g')])
disp(['     which took ' num2str(elapsedcubMCSoboltime) ' seconds and'])
disp(['     ' num2str(paramSobol.n) ' function evaluations to compute'])
disp(' ')


%% Plot the fooling function
figure;
xplot=(0:nplot)/nplot;
fplot=testfun(xplot);
xnode=[(0:0.5:2)*h 2*h+(0:0.25:1)*(1-4*h) 1-(0:0.5:2)*h];
fnode=testfun(xnode);
han=plot(xplot,fplot,'b-',...
    xnode,fnode,'r.',...
    [0 1],exactInteg*[1 1],'k--', ...
    'linewidth',2,'markersize',20);
xlabel('{\it x}')
%ylabel('{\it f(x)}')
legend(han([1 3 2]),{'{\it f(x)}','integral','data'},...
    'location','southwest')
print -depsc FoolQuadFunction.eps


%% Now try the step function
exactInteg=1;
stepfun = @(x,mu,sig,p) mu ...
    + (sig*sqrt((1-p)/p)).*(x<=p) ...
    - (sig*sqrt(p/(1-p))).*(x>p);

p=0.01; 
offset=sqrt(2)/3;
%offset=rand(1);
testfun = @(x) stepfun(mod(x-offset,1),exactInteg,1,p);

%Use quad to find the answer
tic
[quadInteg,quadn]=quad(testfun,0,1,tol);
elapsedquadtime=toc;

%Use quadgk to find the answer
tic
[quadgkInteg,quadgkerr]=quadgk(testfun,0,1,'AbsTol',tol,'RelTol',0);
elapsedquadgktime=toc;

%Use quadl to find the answer
tic
[quadlInteg,quadln]=quadl(testfun,0,1,tol);
elapsedquadltime=toc;

%Use cubMC with iid sampling to find the answer
tic
paramiid.tol=1e-3;
paramiid.sample='iid';
[cubMCiidInteg,paramiid]=cubMC(testfun,[0 1],paramiid);
elapsedcubMCiidtime=paramiid.time;

%Use cubMC with Sobol sampling to find the answer
tic
paramSobol.tol=1e-3;
paramSobol.sample='Sobol';
[cubMCSobolInteg,paramSobol]=cubMC(testfun,[0 1],paramSobol);
elapsedcubMCSoboltime=paramSobol.time;

disp('For the square spike function')
disp(['Exact integral = ' num2str(exactInteg,'%10.15g')])
disp(['Approximate integral using quad = ' num2str(quadInteg,'%10.15g') ...
    ' +/- ' num2str(tol,'%3.1g')])
disp(['     which took ' num2str(elapsedquadtime) ' seconds and'])
disp(['     ' num2str(quadn) ' function evaluations to compute'])
disp(['Approximate integral using quadgk = ' num2str(quadgkInteg,'%10.15g') ...
    ' +/- ' num2str(quadgkerr,'%3.1g')])
disp(['     which took ' num2str(elapsedquadgktime) ' seconds to compute'])
disp(['Approximate integral using quadl = ' num2str(quadlInteg,'%10.15g') ...
    ' +/- ' num2str(tol,'%3.1g')])
disp(['     which took ' num2str(elapsedquadltime) ' seconds and'])
disp(['     ' num2str(quadln) ' function evaluations to compute'])
disp(['Approximate integral using cubMC with iid sampling = '...
    num2str(cubMCiidInteg,'%10.15g') ...
    ' +/- ' num2str(paramiid.tol,'%3.1g')])
disp(['     which took ' num2str(elapsedcubMCiidtime) ' seconds and'])
disp(['     ' num2str(paramiid.n) ' function evaluations to compute'])
disp(['Approximate integral using cubMC with Sobol sampling = '...
    num2str(cubMCSobolInteg,'%10.15g') ...
    ' +/- ' num2str(paramSobol.tol,'%3.1g')])
disp(['     which took ' num2str(elapsedcubMCSoboltime) ' seconds and'])
disp(['     ' num2str(paramSobol.n) ' function evaluations to compute'])
disp(' ')

%% Plot the peak function
figure;
nplot=2000;
xplot=(0:nplot)/nplot;
fplot=testfun(xplot);
han=plot(xplot,fplot,'b-',...
    [0 1],exactInteg*[1 1],'k--', ...
    'linewidth',2,'markersize',20);
xlabel('{\it x}')
%ylabel('{\it f(x)}')
legend(han,{'{\it f(x)}','integral'},'location','best')
axis([0 1 0 12])
print -depsc FoolQuadStepFunction.eps


