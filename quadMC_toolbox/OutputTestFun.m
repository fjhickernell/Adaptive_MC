function OutputTestFun(param)
%   Output the test function chosen by choosetestfun.m

fprintf(1,['For the ' int2str(param.dim) ' dimensional test function named ']);
fprintf(2,[param.fun.funtype ' ']); disp(' ')
if isfield(param,'funDescribe'), fprintf(1,param.funDescribe), end
disp(' ')