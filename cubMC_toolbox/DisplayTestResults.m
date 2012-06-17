%Display TestMCDiffSettings results
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
close all

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
plotTest.nrep=test.nrep;
plotTest.namepref=fun.funtype;
if strcmp(fun.funtype,'step') 
    plotTest.kurtvec=res.exactkurtosis;
    plotTest.namepref=[plotTest.namepref 'd=' int2str(param.dim)];
end



%% Plot iid results
if any(strcmp('iid',test.whichsample))
    plotTest.err=res.iiderr;
    plotTest.time=res.iidtime;
    plotTest.exit=res.iidexit;
    plotTest.kurtmax=res.iidkurtmax;
    plotTest.name=[plotTest.namepref 'iidErrTime'];
    if strcmp(fun.funtype,'step') 
        plotTest.defaultcolor=[0 0 1];
    else
        plotTest.defaultcolor=[1 0 0];
    end
    plotTestcubMC(plotTest,param)
    plotTest=rmfield(plotTest,'kurtmax');
end

%% Plot iid heavy duty results
if any(strcmp('iidheavy',test.whichsample))
    plotTest.err=res.iidheavyerr;
    plotTest.time=res.iidheavytime;
    plotTest.exit=res.iidheavyexit;
    plotTest.kurtmax=res.iidheavykurtmax;
    plotTest.name=[plotTest.namepref 'iidheavyErrTime'];
    if strcmp(fun.funtype,'step') 
        plotTest.defaultcolor=[0 0 1];
    else
        plotTest.defaultcolor=[1 0 0];
    end
    plotTestcubMC(plotTest,param)
    plotTest=rmfield(plotTest,'kurtmax');
end

%% Plot Sobol results
if any(strcmp('Sobol',test.whichsample))
    plotTest.err=res.Sobolerr;
    plotTest.time=res.Soboltime;
    plotTest.exit=res.Sobolexit;
    plotTest.name=[plotTest.namepref 'SobolErrTime'];
    plotTest.defaultcolor=[1 0 0];
    plotTestcubMC(plotTest,param)
end

%% Plot Sobol heavy dutyresults
if any(strcmp('Sobolheavy',test.whichsample))
    plotTest.err=res.Sobolheavyerr;
    plotTest.time=res.Sobolheavytime;
    plotTest.exit=res.Sobolheavyexit;
    plotTest.name=[plotTest.namepref 'SobolheavyErrTime'];
    plotTest.defaultcolor=[1 0 0];
    plotTestcubMC(plotTest,param)
end

%% Plot quad or quadgk results
if any(strcmp('quad',test.whichsample))
    plotTest.err=res.quaderr;
    plotTest.time=res.quadtime;
    plotTest.name=[plotTest.namepref 'quadErrTime'];
    plotTest.defaultcolor=[1 0 0];
    plotTestcubMC(plotTest,param)
end
if any(strcmp('quadgk',test.whichsample))
    plotTest.err=res.quadgkerr;
    plotTest.time=res.quadgktime;
    plotTest.name=[plotTest.namepref 'quadgkErrTime'];
    plotTest.defaultcolor=[1 0 0];
    plotTestcubMC(plotTest,param)
end

