if param.dim==1 && strcmp(param.measure,'uniform');
    %quad answer for test case
    if ~exist('whquad','var'); whquad=''; end
    tic,
    if strcmp(whquad,'gk')
        MatlabCubans=quadgk(@(x) testfun(x'),param.interval(1),param.interval(2),...
            'AbsTol',param.tol); 
    else
        MatlabCubans=quad(@(x) testfun(x'),param.interval(1),param.interval(2),...
            param.tol); 
    end
    MatlabCubtime=toc;
end

if param.dim==2 && strcmp(param.measure,'uniform');
    %dblquad answer for test case
    tic, 
    MatlabCubans=dblquad(@(x,y) testfun([x; repmat(y,1,size(x,2))]'),...
        param.interval(1,1),param.interval(2,1),param.interval(1,2),...
        param.interval(2,2),param.tol); 
    MatlabCubtime=toc;
end

if param.dim==3 && strcmp(param.measure,'uniform');
    %triplequad answer for tes case
    tic, 
    MatlabCubans=triplequad(...
        @(x,y,z) testfun([x; repmat(y,1,size(x,2)); repmat(z,1,size(x,2))]'),...
        param.interval(1,1),param.interval(2,1),param.interval(1,2),...
        param.interval(2,2),param.interval(1,3),param.interval(2,3),...
        param.tol); 
    MatlabCubtime=toc;
end