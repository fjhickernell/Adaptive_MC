disp('Using built-in MATLAB routines')
if param.dim==1
    if strcmp(whquad,'gk')
        fprintf(2,'  quadgk ');
    else
        fprintf(2,'    quad '); 
    end
elseif param.dim==2
    fprintf(2,' dblquad '); 
else
    fprintf(2,'triplelquad '); 
end
disp(['approx. integral = ' num2str(MatlabCubans,'%10.6g')])
if isfield(param,'exactintegral')
    fprintf(2,'            true '); 
    disp(['integral = ' num2str(param.exactintegral,'%10.6g')])
end
fprintf(1,'        Desired '); fprintf(2,'tolerance '); 
disp(['= ' num2str(param.tol)])
if isfield(param,'exactintegral')
    fprintf(1,'             actual '); fprintf(2,'error '); 
    disp(['= ' num2str(abs(param.exactintegral-MatlabCubans),'%10.6g')]); 
end
fprintf(2,'            time '); 
disp(['required = ' num2str(MatlabCubtime) ' seconds'])
disp(' ')