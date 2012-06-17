%Test the cubMC routine

%% Garbage collection and initialization
clear all, close all
format compact

%% Product of exponentials, normal pdf on all reals

% Set up function
fun.funtype='exp';
param.dim=1;
param.measure='uniform';
param.interval=[zeros(1,param.dim);3*ones(1,param.dim)];
%param.measure='normal';
%param.interval=[-inf(1,param.dim);inf(1,param.dim)];
fun.shape=-10./(1:param.dim);
fun.scale=3./(1:param.dim);
fun.addc=1;
fun.overmultc=2;
fun.overaddc=1;
param.impyes=false;
%param.impyes=true;
param.impscale=1.3;
param.impshift=0;
[testfun,param]=choosetestfun(fun,param);
OutputTestFun(param);
OutputTestIntegral(param);

% Evaluate integral
param.tol=1e-3;
param.sample='iid';
[~,param]=cubMC(testfun,param.interval,param);
OutputTestcubMC(param);

param.sample='sobol';
param.scramble=true;
[Q,param]=cubMC(testfun,param.interval,param);
OutputTestcubMC(param);

if param.dim<=3 && strcmp(param.measure,'uniform')
    CubatureMatlab, CubatureMatlabOutput
end