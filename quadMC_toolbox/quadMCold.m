function [Q,param]=quadMC(f,interval,param)
%   QUADMC evaluates a multidimensional integral 
%     by (quasi-)Monte Carlo integration
%     to a specified tolerance (under some assumptions)
%
%          f = the function handle of the integrand
%   interval = the 2 x d multidimensional interval (box) of integration
%      param = a structure with various parameters used in the calculation
%      param.tol      = the error tolerance
%      param.measure  = the name of a measure 
%                         uniform (default) or normal
%                         if normal, the interval is R^d
%      param.sample   = the kind of sampling scheme 
%                         iid (default), lattice, or Sobol
%      param.errmeth  = the kind of error estimation
%                         replications (default) or quasi-standard error
%      param.ndmax    = maximum number of function evaluations
%      param.n0       = the initial sample size
%      param.fudge    = the fudge factor used to estimate the error
%      param.impyes   = perform importance sampling 
%                         (only for normal measure)
%      param.impscale = standard deviation for importance sampling
%      param.impshift = mean shift for importance sampling
%      param.exit     = the state of program when exiting
%      param.fun      = parameters that define the test function (optional)
%
%   The integral to be evaluated takes the form
%       ?   f(x) rho(x) dx
%      R^d
%   where f(x) is the integrand
%         rho(x) is the Radon-Nikodym derivative of the measure
%                = 1 on a finite interval for 'uniform' or
%                = (2*pi)^(-d/2) exp(-x'*x/2) for 'normal'

tstart=tic; %start clock
if nargin<2; interval=[0,1]; end %default interval
if nargin<3; param.tol=1e-2; end %default tolerance
[interval,param]=quadMCparam(interval,param); %check validity of inputs

f=transformIntegrand(f,interval,param); 
%make transformations of the integrand so the sample points don't have to
%   be transformed

if strcmp(param.sample,'iid') %iid sampling
    if strcmp(param.measure,'uniform')
        x=rand(param.n0,param.dim); %uniform random numbers
    else 
        x=randn(param.n0,param.dim); %normal random numbers
    end
    fx=f(x); %evaluate integrand
    sig0=std(fx); %sample standard deviation
    sig0up=param.fudge*sig0; %upper bound on standard deviation
    %estimate sample size needed
    alpha1=1-sqrt(1-param.alpha);
    param.kurtmax=(param.n0-3)/(param.n0-1) ...
        + ((alpha1*param.n0)/(1-alpha1))*(1-1/param.fudge^2)^2;
    A=0.56; %constant in Berry-Esseen inequality
    NBEfun=@(logsqrtn) ...
            normcdf(-exp(logsqrtn).*param.tol/sig0up)...
            +A*param.kurtmax^(3/4).*exp(-logsqrtn)...
            .*(1+exp(logsqrtn).*param.tol/sig0up).^(-3) - alpha1/2;
    sqrtnCLT=norminv(1-alpha1/2)*sig0up/param.tol;
    logsqrtn=fzero(NBEfun,log(sqrtnCLT));
    param.n=ceil(exp(2*logsqrtn));
    %keyboard
    param.n=max(param.n,param.n0);
    if param.n*param.dim>param.ndmax; %too many samples
        param.exit=1;
        [param,Q]=quadMCerr(param,tstart);
        return
    end
    if strcmp(param.measure,'uniform')
        x=rand(param.n,param.dim); %uniform random numbers
    else
        x=randn(param.n,param.dim); %normal random numbers
    end
    fx=f(x); %evaluate integrand
    param.Q=mean(fx); %compute approximate integral as sample mean
    param.n=param.n+param.n0; %total number of samples used

elseif strcmp(param.sample,'Sobol') %Sobol' sampling
    stream=sobolset(param.dim); %create Sobol set
    if param.scramble
        stream=scramble(stream,'MatousekAffineOwen'); %scrambled Sobol set
    end
    stream=qrandstream(stream); %prepare to be used in rand command
    param.n=2^ceil(log2(param.n0)); %initial sample size, a power of 2
    param.npart=2^ceil(log2(param.npart)); %number of pieces of sample
    newn=param.n; %number of samples to take next
    partLength=param.n/param.npart; %number of samples per part
    notDone=true; %whether to keep iterating
    firstTime=true; %whether the iteration is the first time  
    while notDone 
        if param.n*param.dim>param.ndmax; %too many samples needed
            param.exit=1;
            [param,Q]=quadMCerr(param,tstart);
            return
        end
        x=rand(stream,newn,param.dim); %sample next Sobol' points
        fx=f(x); %evaluate integrand
        if firstTime; %partition sampled function values the first time
            fxPart=reshape(fx,partLength,param.npart);
            meanPart=mean(fxPart,1); %sample means of parts
            firstTime=false;
        else %partition sampled function values in successive times
            fxPart=reshape(fx,partLength,param.npart/2);
            meanPart=[(meanPart(1:2:param.npart-1)...
                +meanPart(2:2:param.npart))/2 ...
                mean(fxPart)]; %sample means of parts
        end
        qse=std(meanPart)/sqrt(param.npart); %quasi-standard error for mean
        errEst=param.fudge*qse; %corrected by a fudge factor
        if errEst>param.tol; %not yet satisfied tolerance, don't stop
            newn=param.n; %number of samples to take next
            param.n=2*param.n; %double total sample size
            partLength=partLength*2; %and the length of each part
        else %terminate
            notDone=false;
            param.Q=mean(meanPart); %compute approximate integral
        end
    end
end
param.exit=0; %success!
Q=param.Q; %assign answer
param.time=toc(tstart); %elapsed time
end

function newf=transformIntegrand(oldf,interval,param)
%function to transform the integrand
%    so that the sample points don't have to be transformed

if strcmp(param.measure,'uniform') %uniform measure
    a=interval(1,:); %left endpoint
    b=interval(2,:); %right endpoint
    if all(a~=0) && all(b~=1) %no change needed
        newf=oldf; 
    else %transform points and integrand
        bmina=b-a; %interval width
        volbox=prod(bmina); %volume of the interval
        newf=@(x) oldf(x.*repmat(bmina,size(x,1),1)+repmat(a,size(x,1),1))...
            .*volbox; %stretch and shift, then multiply by volume
    end
elseif strcmp(param.measure,'normal')
    if strcmp(param.sample,'iid') %iid sampling
        if param.impyes
            if all(param.impscale==1) %no expansion or contraction 
                if all(param.impshift==0) %no change needed
                    newf=oldf; 
                else %just a shift
                    newf=@(x) oldf(x+repmat(param.impshift,size(x,1),1)) ...
                        .* exp(-x.*param.impshift'- ...
                        sum(param.impshift.*param.impshift)/2);
                end
            else %expansion or contraction needed
                %keyboard
                if all(param.impshift==0) %no shift
                    newf=@(x) cvfunscaleonly(x,oldf,param.impscale);
                else %need a shift
                    newf=@(x) cvfunscaleshift(x,oldf,...
                        param.impscale,param.impshift);
                end
            end
        else
            newf=oldf;
        end
    else %quasi-Monte Carlo sampling
        if param.impyes
            if all(param.impscale==1) %no expansion or contraction 
                if all(param.impshift==0) %no change needed
                    newf=@(x) oldf(norminv(x)); 
                else %just a shift
                    newf=@(x) cvfunqmcshiftonly(x,oldf,param.impshift);
                end
            else %expansion or contraction needed
                if all(param.impshift==0) %no shift
                    newf=@(x) cvfunscaleonly(norminv(x),oldf,param.impscale);
                else %need a shift
                    newf=@(x) cvfunscaleshift(norminv(x),oldf,...
                        param.impscale,param.impshift);
                end
            end
        else
            newf=@(x) oldf(norminv(x));
        end
    end
end
end

function y=cvfunscaleonly(x,oldf,scale)
%function to perform expansion or contraction only for control variates
    n=size(x,1); %number of points
    %keyboard
    y=oldf(x.*repmat(scale,n,1)) ...
        .* exp(-sum(x.*x.*repmat(scale.*scale-1,n,1),2)/2) ...
        .* prod(scale);
end
        
function y=cvfunscaleshift(x,oldf,scale,shift)
%function to perform expansion or contraction + mean shift
%   for control variates
    n=size(x,1); %number of points
    y=oldf(x.*repmat(scale,n,1)+repmat(shift,n,1)) ...
        .* exp(-((sum(x.*x.*repmat(scale.*scale-1,n,1),2) ...
        + shift*shift'))/2 + x*transpose(shift.*scale)) ...
        .* prod(scale);
end

function y=cvfunqmcshiftonly(x,oldf,shift)
%function to perform mean shift only
%   for control variates with quasi-Monte Carlo sampling
    n=size(x,1); %number of points
    z=norminv(x);  %inverse normal cumulative distribution function
    y=oldf(z+repmat(shift,n,1)) ...
        .* exp(-(shift*shift'/2+z*shift')).* prod(scale).^2;
end
        



    
