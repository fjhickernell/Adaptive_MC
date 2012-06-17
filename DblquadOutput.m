%dblquad answer and output for test cases
tic, 
dbquadans=dblquad(@(x,y) testfun([x; repmat(y,1,size(x,2))]'),...
    param.interval(1,1),param.interval(2,1),param.interval(1,2),param.interval(2,2),param.tol); 
dblquadtime=toc;

%% Print out the approximate integral, error, time, etc.
fprintf(2,' dblquad '); disp(['approx. integral = ' num2str(dbquadans,'%10.6g')])
if isfield(param,'exactintegral')
    fprintf(1,'             actual '); fprintf(2,'error '); disp(['= ' num2str(abs(param.exactintegral-dbquadans),'%10.6g')]); 
end
fprintf(2,'            time '); disp(['required = ' num2str(dblquadtime) ' seconds'])
disp(' ')