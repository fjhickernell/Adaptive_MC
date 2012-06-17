%Run TestcubMC on the GeometricMeanAsianCall
clear all, close all
format compact

fun.funtype='geomean';
fun.S0=100;
fun.K=100;
fun.T=1;
fun.r=0.001;
fun.sigma=0.2;
param.measure='Gaussian';
param.impyes=false;
param.tol=2e-2;
param.n0=1024;
param.dim=8;
param.interval=[-Inf(1,param.dim); Inf(1,param.dim)];
param.sample='Sobol'
[testfun,param]=geomMeanAsianCall(fun,param);

[Q,param]=cubMC(testfun,param.interval,param);
Q
param.time
