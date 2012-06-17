function [testfun,fun,param]=randchoicestep(fun,param,rchparam,irep)

pvec=rchparam.poverall(irep,:);
sigvec=rchparam.sigoverall(irep,:);
fun.shape=pvec;
fun.scale=sigvec./sqrt((1-pvec).*fun.shape);
fun.addc=-sigvec.*sqrt(pvec./(1-pvec));
fun.shift=rand(1,param.dim);
[testfun,param]=choosetestfun(fun,param);
