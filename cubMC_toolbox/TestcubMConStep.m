%Test the new cubMC routine

%% Garbage collection and initialization
clear all, close all
format compact

%% Step function
tstartwhole=tic;

% Set up function
fun.funtype='step';
param.dim=1;
param.measure='uniform';
param.interval=[zeros(1,param.dim); ones(1,param.dim)];
fun.overmultc=1;
fun.overaddc=0;
param.impyes=false;
param.tol=1e-3;

nrep=5;
iidexit=zeros(nrep,1);
iidQ=iidexit;
iiderr=iidexit;
iidtime=iidexit;
iidneval=iidexit;
iidheavyexit=iidexit;
iidheavyQ=iidexit;
iidheavyerr=iidexit;
iidheavytime=iidexit;
iidheavyneval=iidexit;
Sobolexit=iidexit;
SobolQ=iidexit;
Sobolerr=iidexit;
Soboltime=iidexit;
Sobolneval=iidexit;
MatlabQ=iidexit;
Matlaberr=iidexit;
Matlabtime=iidexit;
pmin=1e-5;
pmax=5e-1;
poverall=pmin*(pmax/pmin).^rand(nrep,1);
sigmin=0.1;
sigmax=10;
sigoverall=sigmin*(sigmax/sigmin).^rand(nrep,1);
for irep=1:nrep
    if round(irep/5)==irep/5, irep, end
    fun.shape=poverall(irep);
    fun.scale=sigoverall(irep)/sqrt((1-fun.shape)*fun.shape);
    fun.addc=1-sigoverall(irep)*sqrt(fun.shape/(1-fun.shape));
    fun.shift=rand(1,param.dim);
    [testfun,param]=choosetestfun(fun,param);
    % OutputTestFun(param);
    % OutputTestIntegral(param);

    % Evaluate integral for iid
    param.sample='iid';
    [~,param]=cubMC(testfun,param.interval,param);
    iidkurtmax=param.kurtmax;
    iidexit(irep)=param.exit;
    iidQ(irep)=param.Q;
    iiderr(irep)=abs(param.exactintegral-param.Q);
    iidtime(irep)=param.time;
    iidneval(irep)=param.n;
    %OutputTestcubMC(param);

    % Evaluate integral for heavy duty iid
    nsigold=param.n0;
    param.n0=param.n0*128; %larger n to compute sigma
    param.sample='iid';
    [~,param]=cubMC(testfun,param.interval,param);
    iidheavykurtmax=param.kurtmax;
    iidheavyexit(irep)=param.exit;
    iidheavyQ(irep)=param.Q;
    iidheavyerr(irep)=abs(param.exactintegral-param.Q);
    iidheavytime(irep)=param.time;
    iidheavyneval(irep)=param.n;
    %OutputTestcubMC(param);
    
    % Evaluate integral for Sobol
    param.n0=nsigold;
    param.sample='sobol';
    param.scramble=true;
    [Q,param]=cubMC(testfun,param.interval,param);
    Sobolexit(irep)=param.exit;
    SobolQ(irep)=param.Q;
    Sobolerr(irep)=abs(param.exactintegral-param.Q);
    Soboltime(irep)=param.time;
    Sobolneval(irep)=param.n;
    %OutputTestcubMC(param);

    % Evaluate integral for MATLAB quad
    if param.dim<=3
        CubatureMatlab
        MatlabQ(irep)=MatlabCubans;
        Matlaberr(irep)=abs(param.exactintegral-MatlabCubans);
        Matlabtime(irep)=MatlabCubtime;
        %CubatureMatlabOutput
    end
end

timestamp=datestr(now);
timestamp(timestamp==' ')='_';
timestamp(timestamp==':')='.';
save(['./Results/TestcubMConStepOut' timestamp '.mat']) 
toc(tstartwhole)

DisplayTestStep