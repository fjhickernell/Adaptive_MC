%Test the cubMC routine

%% Garbage collection and initialization
clear all, close all
format compact

%% Product of exponentials, normal pdf on all reals

% Set up function
fun.funtype='step';
param.dim=2;
param.measure='uniform';
param.interval=[zeros(1,param.dim);ones(1,param.dim)];
fun.overmultc=1;
fun.overaddc=1;
param.impyes=false;
myp=0.2*(1:param.dim);
mysig=2;
mykurtvec=1./(myp.*(1-myp))-3
fun.shape=myp;
fun.scale=mysig./sqrt((1-fun.shape).*fun.shape);
fun.addc=-mysig.*sqrt(fun.shape./(1-fun.shape));
fun.shift=0.23;
[testfun,param]=choosetestfun(fun,param);
fun
param.exactvariance
param.exactkurtosis

OutputTestFun(param);
OutputTestIntegral(param);

% Evaluate integral
param.tol=1e-3;
param.n0=1024;
param.sample='iid';
[~,param]=cubMC(testfun,param.interval,param);
OutputTestcubMC(param);

break

param.sample='sobol';
param.scramble=true;
[Q,param]=cubMC(testfun,param.interval,param);
OutputTestcubMC(param);

if param.dim<=3, CubatureMatlab, CubatureMatlabOutput, end