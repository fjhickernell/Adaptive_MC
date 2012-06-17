% Test Geometric Asian Mean
clear all, close all

fun.S0=100;
fun.K=90;
fun.T=1;
fun.d=256;
fun.r=0.03;
fun.sigma=0.5;
param.dim=fun.d;
param.interval=[-Inf(1,param.dim); Inf(1,param.dim)];
param.measure='Gaussian';
fun.funtype='geomean';
[testfun,param]=geomMeanAsianCall(fun,param);
OutputTestFun(param);
OutputTestIntegral(param);

param.tol=1e-1;
param.measure='Gaussian';
param.sample='iid';
[~,param]=cubMC(testfun,param.interval,param);
OutputTestcubMC(param);

param.tol=1e-2;
param.sample='Sobol';
[~,param]=cubMC(testfun,param.interval,param);
OutputTestcubMC(param);
