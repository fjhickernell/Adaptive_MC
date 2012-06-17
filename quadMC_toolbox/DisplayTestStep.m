%Display TestquadMConStep results
close all, clear all
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)

load TestQuadMConStepOut 

plotTest.kurtvec=1./(poverall.*(1-poverall))-3;
plotTest.logerrlo=-5;
plotTest.logerrhi=0;
plotTest.logtimelo=-3;
plotTest.logtimehi=2;
plotTest.errlowlimit=10^plotTest.logerrlo;
plotTest.errhilimit=10^plotTest.logerrhi;
plotTest.timelowlimit=10^plotTest.logtimelo;
plotTest.timehilimit=10^plotTest.logtimehi;
plotTest.linewidth=3;
plotTest.ptsize=500;
plotTest.nrep=nrep;

%% Plot iid results
plotTest.err=iiderr;
plotTest.time=iidtime;
plotTest.exit=iidexit;
plotTest.kurtvec=1./(poverall.*(1-poverall))-3;
plotTest.kurtmax=iidkurtmax;
plotTest.name='iidErrTime';
plotTest.defaultcolor=[0 0 1];
plotTestquadMC(plotTest,param)

%% Plot iid heavy duty results
plotTest.err=iidheavyerr;
plotTest.time=iidheavytime;
plotTest.exit=iidheavyexit;
plotTest.kurtvec=1./(poverall.*(1-poverall))-3;
plotTest.kurtmax=iidheavykurtmax;
plotTest.name='iidheavyErrTime';
plotTest.defaultcolor=[0 0 1];
plotTestquadMC(plotTest,param)

%% Plot Sobol results
plotTest.err=Sobolerr;
plotTest.time=Soboltime;
plotTest.exit=Sobolexit;
plotTest=rmfield(plotTest,'kurtvec');
plotTest=rmfield(plotTest,'kurtmax');
plotTest.name='SobolErrTime';
plotTest.defaultcolor=[1 0 0];
plotTestquadMC(plotTest,param)

%% Plot quad results
plotTest.err=Matlaberr;
plotTest.time=Matlabtime;
plotTest.name='MatlabErrTime';
plotTest.defaultcolor=[1 0 0];
plotTestquadMC(plotTest,param)

