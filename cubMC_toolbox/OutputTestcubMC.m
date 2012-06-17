function OutputTestcubMC(param)
%   Output approximate integral returned by cubMC

fprintf(1,'Using '); fprintf(2,[param.sample ' ']); 
fprintf(1,'sampling');
if param.impyes
    fprintf(1,' with '); fprintf(2,['importance sampling' char(1)]); disp(':')
    disp(['       scale = ' num2str(param.impscale,'%10.6g')]); 
    disp(['       shift = ' num2str(param.impshift,'%10.6g')]); 
else
    disp(':')
end
if param.exit==0;
    fprintf(2,'  cubMC '); disp(['approx. integral = '...
        num2str(param.Q,'%10.6g')])
end
if isfield(param,'exactintegral')
    fprintf(2,'            true '); 
    disp(['integral = ' num2str(param.exactintegral,'%10.6g')])
end
fprintf(1,'        Desired '); fprintf(2,'tolerance '); 
disp(['= ' num2str(param.tol)])
if param.exit==0;
    if isfield(param,'exactintegral')
        fprintf(1,'             actual '); fprintf(2,'error '); 
        disp(['= ' num2str(abs(param.exactintegral-param.Q),'%10.6g')]); 
    end
    fprintf(2,'Number '); disp(['of function values = ' int2str(param.n)])
end
fprintf(2,'            time '); 
disp(['required = ' num2str(param.time) ' seconds, ('...
    num2str(param.time/param.n,'%3.2e') ' seconds/sample)'])
disp(['      for generating samples = ' num2str(param.tsample) ' seconds'])
disp(['    for evaluating integrand = ' num2str(param.tintegrand) ' seconds'])
disp(['        for other operations = ' num2str(param.time ...
    -param.tsample-param.tintegrand) ' seconds'])
disp(' ')