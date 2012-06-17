% Test

p=0.1;
sig=1;
param.tol=1e-6;
fun.shape=p;
fun.scale=sig/sqrt((1-fun.shape)*fun.shape);
fun.addc=1-sig*sqrt(fun.shape/(1-fun.shape));
fun.shift=sqrt(2)/2;
[testfun,param]=choosetestfun(fun,param);
param.sample='iid';
[~,param]=quadMC(testfun,param.interval,param);
param
