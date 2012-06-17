function [funout,param]=importsampGauss(funin,param)
%   Perform a change of variable to do importance sampling for the 
%   Gaussian probability density

param=verifyparam(param,{'impsampscale','impsampshift'},...
    {[1 param.dim],[1 param.dim]},{1,0});
whnotone=param.impsampscale==1;
numnotone=sum(whnotone);
a2min1=param.impsampscale.*param.impsampscale-1;
funout=@(x) funin(repmat(param.impsampscale,1,size(x,1)).*x ...
    + repmat(param.impsampshift,1,size(x,2))) ...
    .* prod(exp(-repmat(a2min1(whnotone),1,numnotone,2)) .* (x(:,whnotone) ...
    + repmat(param.impsampscale(whnotone).*param.impsampshift(whnotone)./a2min1(whnotone),size(x,1)).^2 ...
    - (param.impsampshift(whnotone).^2)/a2min1(whnotone)));