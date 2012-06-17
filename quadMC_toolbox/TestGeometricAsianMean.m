% Test Geometric Asian Mean
clear all, close all

fun.S0=100;
fun.K=100;
fun.T=1;
fun.d=128;
fun.r=0.03;
fun.sigma=0.9;
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
[~,param]=quadMC(testfun,param.interval,param);
OutputTestQuadMC(param);

param.tol=1e-1;
param.sample='Sobol';
[~,param]=quadMC(testfun,param.interval,param);
OutputTestQuadMC(param);
