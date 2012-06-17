function OutputTestIntegral(param)
%   Output the integral defined for evaluation by cubMC.m

fprintf(1,'The integral is taken with respect to the '); 
fprintf(2,[param.measure ' ']);
disp('measure');
fprintf(1,'   on the '); fprintf(2,'interval '); 
numdim=min(6,param.dim);
if param.dim<=numdim
    paramendstring=''; 
else
    paramendstring=' ...'; 
end
disp(['from ' num2str(param.interval(1,1:numdim),'%6.3g') paramendstring]); 
disp(['                     to ' ...
    num2str(param.interval(2,1:numdim),'%6.3g') paramendstring]); 
disp(' ')
