% Example of Monte Carlo
clear, format compact, format long

n=1e6;
d=2;

testfun = @(x,y) exp(-x.^2 + x.*y - y.^2);
tic, exactinteg=dblquad(testfun,0,3,0,3,1e-12), toc

tic, MCappx=9*mean(testfun(3*rand(n,1),3*rand(n,1))), toc

tic, xsobol=scramble(sobolset(d),'MatousekAffineOwen');
QMCappx=9*mean(testfun(3*xsobol(1:n,1),3*xsobol(1:n,2))), toc

errMC=abs(exactinteg-MCappx)/exactinteg
errQMC=abs(exactinteg-QMCappx)/exactinteg