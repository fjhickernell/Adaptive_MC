%Run TestcubMC on the step function
clear all, close all
format compact

fun.funtype='step';
param.dim=1;
param.measure='uniform';
param.interval=[zeros(1,param.dim); ones(1,param.dim)];
fun.overmultc=1;
fun.overaddc=1;
param.impyes=false;
param.tol=1e-3;
param.n0=1024;

test.nrep=1000;
test.howoftenrep=2;
pmin=10^(-5/param.dim);
pmax=5e-1;
test.randch.poverall=pmin*(pmax/pmin).^rand(test.nrep,param.dim);
sigmin=0.1.^(1/param.dim);
sigmax=10.^(1/param.dim);
test.randch.sigoverall=sigmin*(sigmax/sigmin).^rand(test.nrep,param.dim);
test.randchoicefun=@randchoicestep;
test.whichsample={'iid','iidheavy','Sobol','Sobolheavy','quad','quadgk'};
%test.whichsample={'iid','iidheavy','Sobol','Sobolheavy'};
res=TestcubMCDiffSettings(test,fun,param);