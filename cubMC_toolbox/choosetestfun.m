function [testfun,param]=choosetestfun(fun,param)
%   This function chooses and sets up a test function from the parameters
%      input by the user and contained in the structures
%      fun and param
%   fun.funtype         = type of test function
%   param.interval      = domain of test function
%   param.dim           = dimension of the domain
%   param.measure       = measure for integration
%   param.exactintegral = exact value of the integral (scalar)
%   fun.shape           = shape parameter (1 x param.dim)
%   fun.scale           = scale parameter (1 x param.dim)
%   fun.addc            = additive constant (1 x param.dim)
%   fun.overaddc        = overall additive constant (scalar)
%   fun.overmultc       = overall multiplicative constant (scalar)

if nargin < 2; %give the basic default parameters 
    param.interval=[0;1]; %default integration interval
    if nargin < 1; fun.funtype='exp'; end %exponential test function
end
if ~isfield(param,'interval'); param.interval=[0;1]; end %default interval

[~,param]=cubMCparam([],param,'fun'); %check validity of some parameters

if strcmp(fun.funtype,'exp') %exponential test function
    [testfun,param]=makeExpTestFun(fun,param);
elseif strcmp(fun.funtype,'step') %square step test function
    [testfun,param]=makeStepTestFun(fun,param);
else
    error('Function type not recognized')
end

param.fun=fun; %copy of the function parameters
end

%% Exponential Test Function
function [testfun,param]=makeExpTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'shape','scale','addc','overaddc','overmultc'}, ...
    {[1 param.dim],[1 param.dim],[1 param.dim],[1 1],[1 1]}, ...
    {1,1,0,0,1});
testfun=@(x) expfun(x,fun.overaddc,fun.overmultc,fun.addc,fun.scale,...
    fun.shape);

%% Compute exact integral of this function
if strcmp(param.measure,'uniform')
    bmina=param.interval(2,:)-param.interval(1,:);
    param.exactintegral=...
        fun.overaddc.*prod(bmina) ...
        + fun.overmultc.*prod(fun.addc.*bmina ...
        + fun.scale.*(exp(fun.shape.*param.interval(2,:)) ...
        - exp(fun.shape.*param.interval(1,:)))./fun.shape);
else %normal distribution
    param.exactintegral=fun.overaddc ...
        + fun.overmultc.*prod(fun.addc ...
        + fun.scale.*exp(fun.shape.^2/2));
end

%% Prepare text description of the exponential function
if fun.overaddc==0; 
    formulastring='';
else
    formulastring=[num2str(fun.overaddc) ' + ']; 
end
if fun.overmultc~=1; 
    formulastring=[formulastring num2str(fun.overmultc) ' ']; 
end
if param.dim~=1; 
    formulastring=[formulastring 'prod_{j=1}^' int2str(param.dim)];
end
formulastring=[formulastring '['];
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=char(10); 
else
    paramendstring=[' ...' char(10)];
end
if fun.addc==0; 
    paramstring='';
else
    formulastring=[formulastring 'a_j + '];
    paramstring=['    a_j = ' num2str(fun.addc(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.scale~=1; 
    formulastring=[formulastring 'b_j']; 
    paramstring=[paramstring ...
        '    b_j = ' num2str(fun.scale(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.shape==1; 
    formulastring=[formulastring 'exp(x_j)'];
else
    formulastring=[formulastring 'exp(c_j x_j)'];
    paramstring=[paramstring ...
        '    c_j = ' num2str(fun.shape(1:numdim),'%6.3g')...
        paramendstring]; 
end
param.funDescribe=['   f(x) = ' formulastring ']' char(10) paramstring];
end

function f=expfun(x,overaddc,overmultc,addc,scale,shape)
    n=size(x,1);
    f=overaddc+overmultc.*prod(repmat(addc,n,1) ...
        + repmat(scale,n,1).*exp(repmat(shape,n,1).*x),2);
end

%% Step Test Function
function [testfun,param]=makeStepTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'shift','shape','scale','addc',...
    'overaddc','overmultc'}, ...
    {[1 param.dim],[1 param.dim],[1 param.dim],[1 param.dim],...
    [1 1],[1 1]}, ...
    {0.5,0.5,1,1,0,0});
bmina=param.interval(2,:)-param.interval(1,:);
testfun=@(x) stepfun(x,fun.overaddc,fun.overmultc,fun.addc,fun.scale,...
    fun.shape,fun.shift,bmina);

%% Compute exact integral of this function
%   as well as the variance and kurtosis
if strcmp(param.measure,'uniform')
    prodbmina=prod(bmina);
    moment1pc=fun.addc.*bmina+fun.scale.*fun.shape;
    moment2pc=(fun.addc.^2).*bmina+(2.*fun.addc+fun.scale)...
        .*fun.scale.*fun.shape;
    moment3pc=(fun.addc.^3).*bmina+(3.*(fun.addc.^2)...
        +fun.scale.*(3.*fun.addc+fun.scale)).*fun.scale.*fun.shape;
    moment4pc=(fun.addc.^4).*bmina+(4.*(fun.addc.^3)...
        +fun.scale.*(6.*(fun.addc.^2)+fun.scale.*(4.*fun.addc...
        +fun.scale))).*fun.scale.*fun.shape;
    
    moment1=prod(moment1pc);
    moment2=prod(moment2pc);
    moment3=prod(moment3pc);
    moment4=prod(moment4pc);
    param.exactintegral=fun.overaddc.*prodbmina + fun.overmultc.*moment1;
    param.exactvariance=(fun.overmultc.^2).*...
        (moment2-moment1.^2./prodbmina)./prodbmina;
    param.exactkurtosis=(fun.overmultc.^2).*(...
        moment4+moment1.*(-4*moment3+moment1.*(6*moment2...
        -3*moment1.^2)))./(prodbmina.*(param.exactvariance.^2));
    %keyboard
end

%% Prepare text description of the step function
if fun.overaddc==0; 
    formulastring='';
else
    formulastring=[num2str(fun.overaddc) ' + ']; 
end
if fun.overmultc~=1; 
    formulastring=[formulastring num2str(fun.overmultc) ' ']; 
end
if param.dim~=1; 
    formulastring=[formulastring 'prod_{j=1}^' int2str(param.dim)];
end
formulastring=[formulastring '['];
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=char(10); 
else
    paramendstring=[' ...' char(10)];
end
if fun.addc==0; 
    paramstring='';
else
    formulastring=[formulastring 'a_j + '];
    paramstring=['    a_j = ' num2str(fun.addc(1:numdim),'%6.3g')...
        paramendstring]; 
end
if fun.scale~=1; 
    formulastring=[formulastring 'b_j']; 
    paramstring=[paramstring ...
        '    b_j = ' num2str(fun.scale(1:numdim),'%6.3g')...
        paramendstring]; 
end
formulastring=[formulastring ...
    'Indicator_[0,p_j](x - z_j mod(up_j-lo_j)+lo_j'];
paramstring=[paramstring ...
    '    p_j = ' num2str(fun.shape(1:numdim),'%6.3g') paramendstring ...
    '    z_j = ' num2str(fun.shift(1:numdim),'%6.3g') paramendstring ...
    '    lo_j = ' num2str(param.interval(1,1:numdim),'%6.3g') ...
    paramendstring...
    '    hi_j = ' num2str(param.interval(2,1:numdim),'%6.3g') ...
    paramendstring]; 
param.funDescribe=['   f(x) = ' formulastring ']' char(10) paramstring];
end

function f=stepfun(x,overaddc,overmultc,addc,scale,shape,shift,bmina)
    n=size(x,1);
    f=overaddc+overmultc.*prod(repmat(addc,n,1) ...
        + repmat(scale,n,1)...
        .* (mod(x-repmat(shift,n,1),repmat(bmina,n,1))...
        <=repmat(shape,n,1)),2);
end
    

