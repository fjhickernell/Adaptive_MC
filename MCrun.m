function [muhat,ci]=MCrun(testfun,n,mu,sig,sig0,p,wf)
x=rand(n,1);
y=testfun(x,mu,sig,p);
muhat=mean(y);
ci=muhat+(wf*sig0/sqrt(n))*[-1 1];
