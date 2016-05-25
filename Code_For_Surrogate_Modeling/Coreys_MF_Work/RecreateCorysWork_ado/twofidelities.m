function [] = twofidelities()
    gb.x0 = [-2,-2];
    gb.lb = [-2,-2];
    gb.ub= [2,2];
    algorithm = 'sqp';
    highfidelity = @(x)rosenbrock(x); % hfm(x);
    lowfidelity = @(x)sphere(x);%lfm(x);
    
    options.regressionorder = 0;
    options.debug = 0;
    options.hyperinit = [1 1 0.1 0.1];
    optins.hyperspace = [0 0 1 1];
    options.goalfunc = 'mle';
    options.hyperest = 'none';
    options.estvar = 'eigen';
    
    gb.xOpt = gb.x0;
    gb.xCenter = gb.x0;
    gb.x= gb.xOpt;
    [gb.lowOpt,gb.dlowOpt] = lowfidelity(gb.xOpt);
    gb.low= gb.lowOpt;
    gb.dlow = gb.dlowOpt;
    [gb.highOpt,gb.dhighOpt] = highfidelity(gb.xOpt);
    gb.high = gb.highOpt;
    gb.dhigh = gb.dhighOpt;
    wr.wc = [.5 .5];
    tr.size = .25*(gb.ub-gb.lb);
    lowCount = 1;
    highCount = 1;
    converged = 0;
    model.add = [];
    model.mult = [];
    while converged~=1
        [tr,wr,gb] = trmm(tr,wr,gb,model);
        [tr] = adjustmentfactors(tr);
        
        
        if length(tr.add)<=2
            options.regresionorder = 0;
        else
            options.regresionorder = 1;
        end
        model.add = gekPart1(tr.x,tr.add,zeros(size(tr.add)),tr.dadd,zeros(size(tr.dadd)),options);
        
        
        if length(tr.mult)<=2
            options.regresionorder = 0;
        else
            options.regresionorder = 1;
        end
        model.mult = gekPart1(tr.x,tr.mult,zeros(size(tr.mult)),tr.dmult,zeros(size(tr.dmult)),options);
        
        
        [wr] = calculatewc(wr,model);
        [gb.xOpt,~,~,output,~,GRAD] = fmincon(@(x)correctedlf(x,@(x)lowfidelity(x),wr,model),gb.xCenter,[],[],[],[],tr.lb,tr.ub,[],optimset('algorithm',algorithm,'Display','iter','gradobj','on','gradconstr','on'));
        if sum(sum(ismember(gb.x,gb.xOpt),2)==size(gb.xOpt,2))>=1
            gb.xOpt = gb.xOpt - 0.01*sign(GRAD').*gb.xOpt;
            gb.xOpt = min(gb.xOpt,gb.ub);
            gb.xOpt = max(gb.xOpt,gb.lb);
        end
        
        
        n = 25;
        [X,Y] = meshgrid(linspace(tr.lb(1),tr.ub(1),n),linspace(tr.lb(2),tr.ub(2),n));
        CL = zeros(size(X));
        for ii = 1:n
            for jj = 1:n
                CL(ii,jj) = correctedlf([X(ii,jj),Y(ii,jj)],@(x)lowfidelity(x),wr,model);
            end
        end
        surf(X,Y,CL)
        pause(1)
        
        
        lowCount = lowCount+output.funcCount;
        gb.x= [gb.x;gb.xOpt];
        [gb.lowOpt,gb.dlowOpt] = lowfidelity(gb.xOpt);
        gb.low= [gb.low;gb.lowOpt];
        gb.dlow = [gb.dlow;gb.dlowOpt];
        [gb.highOpt,gb.dhighOpt] = highfidelity(gb.xOpt);
        gb.high = [gb.high;gb.highOpt];
        gb.dhigh = [gb.dhigh;gb.dhighOpt];
        lowCount = lowCount+1;
        highCount = highCount+1;
        epsilon = abs((gb.highOpt-gb.highCenter)/gb.highOpt);
        if epsilon<=0.001&&epsilonOld<0.001&&tr.accept==1
            converged = 1;
        end
        epsilonOld = epsilon;
    end
    gb.xOpt
    gb.highOpt
    lowCount
    highCount
end

function [f,df] = hfm(x)
        f = ((2*x-2).^2.*sin(12*x-4)).*(5*x.^3-3*x.^2+7*x+1)+(-4*x.^2+2*x-1);
    if nargout == 2
        df = sin(12*x-4).*(2*x-2).^2.*(15*x.^2-6*x+7)-8*x+12*cos(12*x-4).*(2*x-2).^2.*(5*x.^3-3*x.^2+7*x+1)+sin(12*x-4).*(8*x-8).*(5*x.^3-3*x.^2+7*x+1)+2;
    end
end

function [f,df] = mhfm(x)
        f = ((2*x-2).^2.*sin(12*x-4)).*(5*x.^3-3*x.^2+7*x+1);
    if nargout == 2
        df = sin(12*x-4).*(2*x-2).^2.*(15*x.^2-6*x+7)+12*cos(12*x-4).*(2*x-2).^2.*(5*x.^3-3*x.^2+7*x+1)+sin(12*x-4).*(8*x-8).*(5*x.^3-3*x.^2+7*x+1);
    end
end

function [f,df] = mlfm(x)
        f = ((2*x-2).^2.*sin(12*x-4))+(-4*x.^2+2*x-1);
    if nargout == 2
        df = sin(12*x-4).*(8*x-8)-8*x+12*cos(12*x-4).*(2*x-2).^2+2;
    end
end

function [f,df] = lfm(x)
        f = (2*x-2).^2.*sin(12*x-4);
    if nargout == 2
        df = sin(12*x-4).*(8*x-8)+12*cos(12*x-4).*(2*x-2).^2;
    end
end

function [f,df] = rosenbrock(x)
    f = sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1:end-1)-1).^2);
    if nargout > 1
        n = length(x);
        df(1,2:n-1) = -400*reshape(x(2:n-1),n-2,1).*(reshape(x(3:n),n-2,1)...
            -reshape(x(2:n-1),n-2,1).^2)+2*(reshape(x(2:n-1),n-2,1)-1)+...
            200*(reshape(x(2:n-1),n-2,1)-reshape(x(1:n-2),n-2,1).^2);
        df(1,1) = -400*(x(2)-x(1)^2)*x(1)+2*(x(1)-1);
        df(1,n) = 200*(x(n)-x(n-1)^2);
    end
end

function [f,df,H] = branin(x)
    a = 1;
    b = 5.1/4/pi^2;
    c = 5/pi;
    r = 6;
    s= 10;
    t = 1/8/pi;

    f = s + a*(b*x(:,1).^2 - c*x(:,1) + r - x(:,2)).^2 - s*cos(x(:,1))*(t - 1);
    if nargout > 1
        df(:,1) = s*sin(x(:,1))*(t - 1) - 2*a*(c - 2*b*x(:,1)).*(b*x(:,1).^2 - c*x(:,1) + r - x(:,2));
        df(:,2) = -a*(2*b*x(:,1).^2 - 2*c*x(:,1) + 2*r - 2*x(:,2));
    end
end

function [f,df,H] = sphere(x)
    f = sum(100*x.^2);
    if nargout > 1
        df = 200*x;
    end
end

function [f,df] = correctedlf(x,lf,wr,model)
    if isa(lf,'function_handle')
        if nargout == 1
            low = lf(x);
        else
            [low,dlow] = lf(x);
        end
    else
        low = lf(1);
        dlow = lf(2:end);
    end
    
    if nargout == 1
        beta = additivecorrection(x,model.add);
        gamma = multiplicativecorrection(x,model.mult);
    else
        [beta,dbeta] = additivecorrection(x,model.add);
        [gamma,dgamma] = multiplicativecorrection(x,model.mult);
    end
    
    f = wr.wc(1)*(low+beta)+wr.wc(2)*(low*gamma);
    if nargout == 2
        df = wr.wc(1)*(dlow+dbeta)+wr.wc(2)*(low*dgamma+gamma*dlow);
    end
end

function [beta,dbeta] = additivecorrection(x,model)
    if nargout == 1
        [beta,~,~,~,~] = gekPart2(model,x);
    else
        [beta,~,~,~,~] = gekPart2(model,x);
        dbeta = zeros(size(x));
        for ii = 1:length(x)
            tmpX = x;
            tmpX(ii) = tmpX(ii)+.0001;
            [betaFor,~,~,~,~] = gekPart2(model,tmpX);
            dbeta(ii) = (betaFor-beta)/.0001;
        end
    end
end

function [gamma,dgamma] = multiplicativecorrection(x,model)
    if nargout == 1
        [gamma,~,~,~,~] = gekPart2(model,x);
    else
        [gamma,~,~,~,~] = gekPart2(model,x);
        dgamma = zeros(size(x));
        for ii = 1:length(x)
            tmpX = x;
            tmpX(ii) = tmpX(ii)+.0001;
            [gammaFor,~,~,~,~] = gekPart2(model,tmpX);
            dgamma(ii) = (gammaFor-gamma)/.0001;
        end
    end
end

function [tr,wr,gb] = trmm(tr,wr,gb,model)
    r = [0,1e-4,.8];
    s = [.5,.75,1.25,2.5];
    if size(gb.x,1)==1
        rho = 1;
    else
        rho = (gb.highCenter-gb.highOpt)/(gb.highCenter-correctedlf(gb.xOpt,[gb.lowOpt,gb.dlowOpt],wr,model));
    end
    if rho>r(1)
        gb.xCenter = gb.xOpt;
        gb.highCenter = gb.highOpt;
        fprintf('Accept Step - ')
        tr.accept = 1;
    else
        fprintf('Reject Step - ')
        tr.accept = 0;
    end
    if rho<=r(1)
        scale = s(1);
        fprintf('Inaccurate\n\n')
    elseif rho<=r(2)
        scale = s(2);
        fprintf('Marginally Accurate\n\n')
    elseif rho<=r(3)
        scale = s(3);
        fprintf('Moderately Accurate\n\n')
    else
        scale = s(4);
        fprintf('Accurate\n\n')
    end
    tr.size = scale*tr.size;
    tr.lb = max(gb.lb,gb.xCenter-tr.size);
    tr.ub = min(gb.ub,gb.xCenter+tr.size);
    
    indLB = ge(gb.x,ones(size(gb.x,1),1)*tr.lb);
    indUB = le(gb.x,ones(size(gb.x,1),1)*tr.ub);
    ind = intersect(find(sum(indLB,2)==size(gb.x,2)),find(sum(indUB,2)==size(gb.x,2)));

    tr.x = gb.x(ind,:);
    tr.low = gb.low(ind);
    tr.dlow = gb.dlow(ind,:);
    tr.high = gb.high(ind);
    tr.dhigh = gb.dhigh(ind,:);
    
    wr.lb = max(gb.lb,gb.xCenter-1.5*tr.size);
    wr.ub = min(gb.ub,gb.xCenter+1.5*tr.size);
    
    indLB = ge(gb.x,ones(size(gb.x,1),1)*wr.lb);
    indUB = le(gb.x,ones(size(gb.x,1),1)*wr.ub);
    ind = intersect(find(sum(indLB,2)==size(gb.x,2)),find(sum(indUB,2)==size(gb.x,2)));

    wr.x = gb.x(ind,:);
    wr.low = gb.low(ind);
    wr.dlow = gb.dlow(ind,:);
    wr.high = gb.high(ind);
    wr.dhigh = gb.dhigh(ind,:);

    indLB = le(wr.x,ones(size(wr.x,1),1)*tr.lb);
    indUB = ge(wr.x,ones(size(wr.x,1),1)*tr.ub);
    ind = intersect(find(sum(indLB,2)==size(wr.x,2)),find(sum(indUB,2)==size(wr.x,2)));

    wr.x = wr.x(ind,:);
    wr.low = wr.low(ind);
    wr.dlow = wr.dlow(ind,:);
    wr.high = wr.high(ind);
    wr.dhigh = wr.dhigh(ind,:);
end

function [tr] = adjustmentfactors(tr)
    tr.add = tr.high-tr.low;
    tr.dadd = tr.dhigh-tr.dlow;
    tr.mult = tr.high./tr.low;
    tr.dmult = (tr.low*ones(1,size(tr.x,2)).*tr.dhigh-tr.high*ones(1,size(tr.x,2)).*tr.dlow)./(tr.low.^2*ones(1,size(tr.x,2)));
    idx = tr.low==0;
    tr.mult(idx) = [];
    tr.dmult(idx,:) = [];   
end

function [wr] = calculatewc(wr,model)
    n = size(wr.x,1);
    if n>=1
        for ii = 1:n
            ADD(ii,1) = correctedlf(wr.x(ii,:),[wr.low(ii),wr.dlow(ii,:)],[1 0],model);
            MULT(ii,1) = correctedlf(wr.x(ii,:),[wr.low(ii),wr.dlow(ii,:)],[0 1],model);
        end
        errorADD = wr.high-ADD;
        errorMULT = wr.high-MULT;
        sig2ADD = sum(errorADD.^2)/n;
        sig2MULT = sum(errorMULT.^2)/n;
        psiADD = (1/(2*pi*sig2ADD))^(n/2)*exp(-n/2);
        psiMULT = (1/(2*pi*sig2MULT))^(n/2)*exp(-n/2);
        psi = [psiADD;psiMULT];
        wr.wc = [weightCoeffs(1)*psi(1), weightCoeffs(2)*psi(2)]/(weightCoeffs*psi);
        if (isnan(wr.wc(1))||isinf(wr.wc(1)))&&~(isnan(wr.wc(2))||isinf(wr.wc(2)))
            wr.wc = [.00001 .99999];
        elseif (isnan(wr.wc(2))||isinf(wr.wc(2)))&&~(isnan(wr.wc(1))||isinf(wr.wc(1)))
            wr.wc = [.99999 .00001];
        elseif (isnan(wr.wc(2))||isinf(wr.wc(2)))&&(isnan(wr.wc(1))||isinf(wr.wc(1)))
            wr.wc = [.5 .5];
        end
    else
        wr.wc = [.5 .5];
    end
end

function gekmodel = gekPart1(xi,y,erry,grad,errgrad,options)
%% GEKPART1 Gradient-Enhanced Kriging - Part 1
% 
% gekmodel = ...
%    gekPart1(xi,y,erry,grad,errgrad,options)
%
% GEK can be computationally expensive. In many cases, the data is 
% available before xiout has been determined. In these cases, one might 
% want to use gekPart1 and gekPart2.
%
% gekPart1 runs the expensive computations, such as hyperparameter
% estimation and computation of the initial state.
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% y - observed QoI, size N x 1
% erry - measurement error in y, size N x 1
% grad - observed gradients, size N x ndim
% errgrad - measurement error in grad, size N x ndim
% xiout - location of predictions, size n x ndim
% options - GEK options, refer to defaultopts
%
% NOTE: specify '' for grad and errgrad if no gradient information is 
% available
%
% == OUTPUT ==
% gekmodel - struct containing model info
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

%% === add default options ===
options = ...
    defaultopts(options,xi);

if options.debug > 0
    fprintf('\n=== gek part 1 ===\n')
end

%% === normalize ===
[yn errn norming N normreport] = ...
    normalize(xi,y,erry,grad,errgrad,options);

%% === estimate hyperparameters ===
% this is still incomplete
[hypers hypreport] = ...
    hyperestimate(xi,yn,errn,options);

%% === find initial state ===
[state0 A initreport] = ...
    initialstate(xi,yn,errn,hypers,options);

gekmodel.xi = xi;
gekmodel.state0 = state0;
gekmodel.hypers = hypers;
gekmodel.A = A;
gekmodel.options = options;
gekmodel.N = N;
gekmodel.norming = norming;
end

function [xout varxout gradout vargradout report] = gekPart2(gekmodel,xiout)
%% GEKPART2 Gradient-Enhanced Kriging - Part 2
% 
% [xout varxout gradout vargradout report] = ...
%    gekPart2(gekmodel,xiout)
%
% GEK can be computationally expensive. In many cases, the data is 
% available before xiout has been determined. In these cases, one might 
% want to use gekPart1 and gekPart2.
%
% gekPart2 only runs the (less expensive) compations which require xiout.
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% gekmodel - struct containing model info
% xiout - location of predictions, size n x ndim
%
% NOTE: specify '' for grad and errgrad if no gradient information is 
% available
%
% == OUTPUT ==
% xout - predicted mean of QoI, size n x 1
% varxout - predicted variance of QoI, size n x 1
% gradout - predicted gradients, size n x ndim
% vargradout - predicted variance of gradients, size n x ndim
% report - << under construction >>
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

%%
if gekmodel.options.debug > 0
    fprintf('\n=== gek part 2 ===\n')
end

%% === make prediction ===
[xoutn varxoutn predreport] = ...
    prediction(gekmodel.xi,gekmodel.state0,gekmodel.hypers,...
        gekmodel.A,xiout,gekmodel.options);

%% === denormalize prediction ===
[xout varxout gradout vargradout denormreport] = ...
    denormalize(xiout,xoutn,...
        varxoutn,gekmodel.norming,gekmodel.options);

%% === create final report ===
report = 'to do';
%    createreport(normreport,hypreport,initreport,predreport,denormreport);
end

function options = defaultopts(options,xi)
%% DEFAULTOPTS Complements options for GEK
% 
% options = defaultopts(options,xi)
%
% Complements options with default GEK options. Please refer to code below
% for more details.
%
% == INPUT ==
% options - struct with options
% xi - location of observations, size N x ndim
%
% == OUTPUT ==
% options - struct with options
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

[n ndim] = size(xi);

%% general
if isfield(options,'debug') == 0
    options.debug = 0; % 0 for no debug
                       % 1 for command line feedback
                       % 2 for cl feedback & plots
                       % -1 for test-phase code
end


%% normalize.m
if isfield(options,'regressionorder') == 0
    options.regressionorder = 0;
end

if isfield(options,'nregresscoeff') == 0
    options.nregresscoeff = 2 + ndim*options.regressionorder;
end

if isfield(options,'stdlowerbound') == 0
    options.stdlowerbound = 1e-3;
end

%% hyperestimate.m
if isfield(options,'hyperinit') == 0
    rangeinit = 0.3 * ( max(xi) - min(xi) );
    options.hyperinit = [1 1 rangeinit];
end

if isfield(options,'hyperspace') == 0
    options.hyperspace = [0 0 ones(1,ndim)];
end

if isfield(options,'goalfunc') == 0
    options.goalfunc = 'mle';
end

if isfield(options,'hyperest') == 0
    options.hyperest = 'fmin';
end

if isfield(options,'fminopts') == 0
    options.fminopts = ...
        optimset('Display','off','TolFun',1e-20,...
            'TolX',1e-3,'MaxIter',1e3);
end

if isfield(options,'brutesize') == 0
    options.brutesize = 256;
end

if isfield(options,'brutefactorlimits') == 0
    options.brutefactor = sqrt(10);
end

%% predict.m
if isfield(options,'estvar') == 0
    options.estvar = 'yes';
end

if isfield(options,'predicttype') == 0
    options.predicttype = 'predict';
end

if isfield(options,'maxArraySize') == 0
    %if exist('memory') == 0
        options.maxArraySize = 1e8;
    %else
    %    user = memory;
    %    safety = 1.3;
    %    options.maxArraySize = user.MaxPossibleArrayBytes / 8 / safety;
    %end
end
end

function [yn errn norming N normreport] = normalize(xi,y,erry,grad,errgrad,options)
%% NORMALIZE Normalize data for GEK
% 
% [yn errn norming N normreport] = ...
%    normalize(xi,y,erry,grad,errgrad,options)
%
% Normalizes data for GEK. In case of gradients, compiles value and
% gradient information to single data vector.
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% y - observed QoI, size N x 1
% erry - measurement error in y, size N x 1
% grad - observed gradients, size N x ndim
% errgrad - measurement error in grad, size N x ndim
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% yn - normalized data vector, size N x (1+ndim)
% errn - normalized error vector, size N x (1+ndim)
% norming - struct containing normalization info
% normreport - << under construction >>
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

tic
%% == general ==
[ny ndim] = size(xi);
N = ny;

if options.debug > 0
    fprintf('gek -> normalizing ...\n')
end

%% == normalize data ==
if options.regressionorder == 0
    norming.c = [mean(y) ; zeros(ndim,1)]; % coefficients for linear drift function
elseif options.regressionorder == 1
    norming.c = robustfit(xi,y,'ols',''); % coefficients for linear drift function
else
    warning('higher order regression not supported')
end
%norming.c;
ylin = norming.c(1) + xi*norming.c(2:end);
gradlin = ones(ny,1)*norming.c(2:end)';
stdev = std(y-ylin);
if stdev > options.stdlowerbound
    norming.st = stdev;
else
    norming.st = options.stdlowerbound;
    if options.debug > 0, warning('standard deviation reached lower bound'), end
end

% normalization
ynn = (y-ylin)/norming.st;
errynn = erry/norming.st;
if ~strcmp(grad,'')
    gradnn = (grad-gradlin)/norming.st;
    errgradnn = errgrad/norming.st;
    N = (ndim+1)*N;
end

% compile data
if strcmp(grad,'')
    yn = ynn;
    errn = errynn;
else
    gradnn = reshape(gradnn,ny*ndim,1);
    errgradnn = reshape(errgradnn,ny*ndim,1);
    yn = [ynn ; gradnn];
    errn = [errynn ; errgradnn];
end

%% == report ==
normreport.cpu = toc;
end

function [state0 A initreport] = initialstate(xi,yn,errn,hypers,options)
%% INITIALSTATE Initial state for GEK
% 
% [state0 A initreport] = ...
%    initialstate(xi,yn,errn,hypers,options)
%
% Computes the initial state: 
% 
%    state0 = A\yn, 
%
% where A = R + H'PH. Note that the initial state can contain value and
% gradient information, since yn is a (compiled) normalized data vector
%
% Calculation of the initial state is a expensive operation (deblurr or 
% pull-back). The initial state can later be used for cheap prediction
% of the QoI at xiout (blurr or push-forward or diffusion)
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% yn - normalized compiled data vector, size N x (1+ndim)
% errn - normalized error vector, size N x (1+ndim)
% hypers - optimized hyperparameters, size 1 x (ndim + 2)
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% state0 - initial state, size N x (1+ndim)
% A - gain matrix, size N(1+ndim) x N(1+ndim)
% initreport - << under construction >>
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

%% == check for covariables ==
tic
[n1 m1] = size(xi);
[n2 m2] = size(yn);
if n1 == n2
    covars = 'no';
else
    covars = 'yes';
end

if options.debug > 0
    fprintf('gek -> finding initial state ...\n')
end

%% == split hypers ==
errnf = nan(n2,1);
if n1 == n2
    errnf = hypers(1) * errn;
else
    errnf(1:n1) = hypers(1) * errn(1:n1);
    errnf(n1+1:end) = hypers(2) * errn(n1+1:end);
end

range = hypers(3:end);

%% == determine initial state ==
   
P = corrmatrix(xi,xi,range,'analyze',covars,options);             % prior
R = diag(errnf.^2);                   % likelihood
A = R + P;                           % gain (is output for later use)
state0 = A\yn;                       % cholesky -> backward substitution

%% == check the accuracy of the initial state ==
prec = eps;
accur = std(A*state0-yn);
if options.debug > 0
    fprintf(['   residual : ' num2str(accur) '\n'])
    if accur > sqrt(prec);
        warning(['the initial state for size(xi) = [ ' num2str(size(xi)) ' ] is quite fuzzy'])
    end
end

%% == report ==
initreport.cpu = toc;
initreport.accuracy = accur;

if options.debug > 2
    if strcmp(covars,'no'), figure(11), else figure(11), end
    subplot(2,1,1),contourf(A,100,'edgecolor','none'), colorbar
    if strcmp(covars,'no'), title('kriging A'), else title('gek A'), end
    subplot(2,1,2),contourf(chol(A),100,'edgecolor','none'), colorbar
    if strcmp(covars,'no'), title('kriging chol(A)'), else title('gek chol(A)'), end
end
end

function P = corrmatrix(xi1,xi2,hypers,type,covars,options)
%% CORRMATRIX Corrmatrix (for GEK)
%
% P = corrmatrix(xi1,xi2,hypers,type,covars,options)
%
% Correlation matrix for Gradient-Enhanced Kriging
%
% see also HYPERESTIMATE, INITIALSTATE, PREDICT

%% == general ==
[nxi1 ndim] = size(xi1);
[nxi2 ndim] = size(xi2);
                
%% == lag ==

% check memory and preallocate lag
if options.debug > 0
    mag = nxi1*nxi2*ndim; % required memory
    if mag > options.maxArraySize
        warning('memory might be insufficient to construct lag matrix')
    end
end
H = zeros(nxi1,nxi2,ndim); Hn = H; % preallocate

% slow implementation
% tic
% for i = 1:nxi1
%     for j = 1:nxi2
%         for k = 1:ndim
%             H(i,j,k) = xi1(i,k) - xi2(j,k);
%             Hn(i,j,k) = ( xi1(i,k) - xi2(j,k) ) / hypers(k) ;
%         end
%     end
% end
% toc

% fast implementation
%tic
for k = 1:ndim
    XI1 = xi1(:,k)*ones(1,nxi2);
    XI2 = ones(nxi1,1)*xi2(:,k)';
    H(:,:,k) = XI1 - XI2;
    Hn(:,:,k) = (XI1 - XI2)/hypers(k);
end
%toc

%% == correlation matrix ==

% determine number of rows of matrix blocks
% for the analysis we need the full matrix
% for the prediction we only need the top row
if strcmp(type,'analyze') && strcmp(covars,'yes')
    rows = ndim + 1; % construct full matrix: all rows of blocks
elseif strcmp(type,'predictplus')
    rows = ndim + 1;
else
    rows = 1; % construct only top row of blocks
end

if strcmp(covars,'yes')
    cols = ndim + 1;
elseif strcmp(covars,'no')
    cols = 1;
end

% check memory and preallocate correlation matrix
if options.debug > 0
    mag = nxi1*rows*nxi2*(ndim+1); % required memory
    if mag > options.maxArraySize
        warning('memory might be insufficient to construct correlation matrix')
    end
end
P = zeros(nxi1*rows,nxi2*cols);

% precompute matrices
Hn2 = sum(Hn.^2,3);
EXP = exp(-0.5*Hn2);

% fill blocks
for i = 1:rows
    for j = 1:cols
        if i == 1 && j == 1 % no derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                EXP;
        elseif i == 1 % first derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                H(:,:,j-1) .* hypers(j-1)^-2 .* EXP;
        elseif j == 1 % first derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                - H(:,:,i-1) .* hypers(i-1)^-2 .* EXP;
        elseif i == j % second derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                hypers(i-1)^-2 * EXP - ...
                H(:,:,i-1).^2 .* hypers(i-1)^-4 .* EXP;
        else % two partial derivatives
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                - H(:,:,i-1) .* hypers(i-1)^-2 .* ...
                H(:,:,j-1) .* hypers(j-1)^-2 .* ...
                EXP;
        end
    end
end
end

function measure = goalfunction(loghyper,hypid,hypinit,xi,yn,errn,type,gradients)
%% GOALFUNCTION Goalfunction (for GEK)
%
% measure = goalfunction(loghyper,hypid,hypinit,xi,yn,errn,type,gradients)
%
% Goal function for optimization of hyperparameters. Calculates negative
% log likelihood.
%
% see also HYPERESTIMATE

mleopts.debug = 0;

N = size(xi,1);
N2 = size(yn,1);

%% reconstruct hyper
fullhyper = hypinit;
fullhyper(hypid) = exp(loghyper);

if N == N2
    erryfactor = fullhyper(1);
else
    erryfactor = fullhyper(1);
    errgradfactor = fullhyper(2);
end
range = fullhyper(3:end);

errnf = nan(N2,1);
if N == N2
    errnf = erryfactor * errn;
else
    errnf(1:N) = erryfactor * errn(1:N);
    errnf(N+1:end) = errgradfactor * errn(N+1:end);
end

P = corrmatrix(xi,xi,range,'analyze',gradients,mleopts);             % prior
R = diag(errnf.^2);                       
A = R + P;

switch type
    case 'mle'
        measure = sum(log(eig(A))) + yn'*(A\yn);
    case 'mle_efficient'
        % more efficient, less robust
        [U p] = chol(A);
        if 0 == 0
           % this is more efficient
           measure = sum(log(eig(U))) + yn'*(U'\(U\yn));
        else
           % this is more robust
           measure = sum(log(eig(A))) + yn'*(A\yn);
        end
    case 'mle_selected'
        % ???? this might improve results for small N ????
        f = ones(size(yn,1),1);
        if size(yn,1) > size(xi,1)
            f(size(xi,1)+1:end) = 0; % do we want t
        end
        measure = sum(log(eig(A))) + (f.*yn)'*(A\yn);
    otherwise
        warning('this type of goal function is not supported')
end
end

function [hypers hypreport] = hyperestimate(xi,yn,errn,options)
%% HYPERESTIMATE Estimate hyperparameters for GEK
% 
% [hypers hypreport] = ...
%    hyperestimate(xi,yn,errn,options)
%
% Optimizes the hyperparameters for GEK, using a Maximum Likelihood 
% Estimate (MLE). 
%
% Note that in the code there are ndim + 2 hyperparameters, i.e. one
% erryfactor, one errgradfactor, and ndim correlation ranges. By default,
% GEK does not optimize for the err factors, and optimizes for all
% correlation ranges.
%
% Optimizing the hyperparameters is often the most expensive part of GEK in
% terms of CPU time, so it might pay to provide good initiall guesses and
% consider the relevant options, such as:
%
% options.debug = 0; prints optimization results <default>
% options.debug = 1; prints optimization results
% options.debug = 2; prints & plots optimization results
%
% options.hyperest = 'fmin'; uses fminsearch optimization <default>
% options.hyperest = 'brute'; uses brute force optimization
%
% options.hyperinit = [1 1 0.3 0.4]; this is an example setting for the
%  initial guess of the hyperparameters for a 2-d case, where erryfact = 1, 
%  errgradfact = 1, correlation range1 = 0.3, and range2 = 0.4
%  <default setting depends on N and on domain size>
%
% options.hyperspace = [0 0 1 1]; 2-d example, optimizes only for both 
%  correlation ranges <default>
% options.hyperspace = [0 1 0 1]; 2-example, optimizes only for the
%  errgradfact and for correlation range2
%
% options.brutesize = 256; # samples in brute optimization <default>
% 
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% yn - normalized data vector, size N x (1+ndim)
% errn - normalized error vector, size N x (1+ndim)
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% hypers - optimized hyperparameters, size 1 x (ndim + 2)
% hypreport - << under construction >>
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

%% == general ==
[ny ndim] = size(xi);
[nyn dum] = size(yn);

if nyn > ny
    gradients = 'yes';
else
    gradients = 'no';
end

tic

%% == estimate hypers ==
switch options.hyperest
    case 'none'
        if options.debug > 0
            fprintf('gek -> fixed hyperparameters \n')
        end
        hypers = options.hyperinit;
        erFact = 'empty';
    case 'brute'
        if options.debug > 0
            fprintf('gek -> estimating hyperparameters (brute force) ... \n')
        end
        hypspace = options.hyperspace;
        hypinit = options.hyperinit;
        hypid = find(hypspace~=0);
        hypdim = length(hypid);
        hyp0 = lhsdesign(options.brutesize,hypdim,'iter',256);
        hyp0 = (options.brutefactor).^(2*hyp0-1);
        hyp = repmat(hypinit(hypid),options.brutesize,1) .* hyp0;
        loghyp = log(hyp);
        for i = 1:options.brutesize
            negloglike(i) = ...
                goalfunction(loghyp(i,:),hypid,hypinit,...
                    xi,yn,errn,options.goalfunc,gradients);
        end
        [m im] = min(negloglike);
        hypers = hypinit;
        hypers(hypid) = exp(loghyp(im,:));
        if hypdim == 1 && options.debug > 1
            figure(10)
            loglog(exp(loghyp),exp(-negloglike),'b.')
            hold on 
            loglog(exp(loghyp(im)),exp(-m),'ko')
            hold off
            xlabel(['hyper ' num2str(hypid)])
            title('likelihood')
        elseif hypdim == 2 && options.debug > 1
            figure(10)
            scatter(log10(exp(loghyp(:,1))),log10(exp(loghyp(:,2))),...
                3,exp(-negloglike))
            hold on
            scatter(log10(exp(loghyp(im,1))),log10(exp(loghyp(im,2))),...
                40,'ko','filled')
            hold off
            xlabel(['log10 hyper ' num2str(hypid(1))])
            ylabel(['log10 hyper ' num2str(hypid(2))])
            colorbar
            title('likelihood')
        end
    case 'fmin'
        hypspace = options.hyperspace;
        hypinit = options.hyperinit;
        hypid = find(hypspace~=0);
        [loghypfmin,fval,exitflag,output] = ...
            fminsearch(@(loghyp)goalfunction(loghyp,hypid,hypinit,...
                    xi,yn,errn,options.goalfunc,gradients),...
                        log(hypinit(hypid)),options.fminopts);
        hypers = hypinit;
        hypers(hypid) = exp(loghypfmin);
    otherwise
        error('this type of minimization of the goalfunction is not supported')
end

if options.debug > 0
    fprintf(['     erryfact    : ' num2str(hypers(1)) '\n'])
    fprintf(['     errgradfact : ' num2str(hypers(2)) '\n'])
    fprintf(['     range (the) : ' num2str(hypers(3:end)) '\n'])
end

%% == report ==
hypreport.cpu = toc;
end

function [xout varxout gradout vargradout denormreport] = denormalize(xiout,xoutn,varxoutn,norming,options)
%% DENORMALIZE De-normalize data for GEK
% 
% [xout varxout gradout vargradout denormreport] = ...
%    denormalize(xiout,xoutn,varxoutn,norming,options)
%
% De-normalizes data for GEK. In case of gradients, de-compiles value and
% gradient information individual vectors.
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xiout - location of predictions, size n x ndim
% xoutn - predicted normalized mean, size n x 1
% varxoutn - predicted normalized variance, size n x 1
% norming - struct containing normalization info
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% xout - predicted mean of QoI, size n x 1
% varxout - predicted variance of QoI, size n x 1
% gradout - predicted gradients, size n x ndim
% vargradout - predicted variance of gradients, size n x ndim
% denormreport - << under construction >>
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

tic;
%% == general ==
[nx ndim] = size(xiout);

if options.debug > 0
    fprintf('gek -> denormalizing ...\n')
end

%% == denormalize ==
% mean
ylin = norming.c(1) + xiout*norming.c(2:end);
xout = norming.st*xoutn(1:nx) + ylin;

% variance
if strcmp(options.estvar,'yes')
    varxout = norming.st^2 * varxoutn(1:nx);
elseif strcmp(options.estvar,'cheap')
    varxout = norming.st^2 * varxoutn(1);
else
    varxout = '';
end

% gradients
if strcmp(options.predicttype,'predictplus')
    gradout = norming.st*xoutn(nx+1:end);
    gradout = reshape(gradout,nx,ndim);
    vargradout = norming.st*varxoutn(nx+1:end);
%     vargradout = reshape(vargradout,nx,ndim);
else
    gradout = NaN;
    vargradout = NaN;
end

%% == report ==
denormreport.cpu = toc;
end

function [xoutn varxoutn predreport] = prediction(xi,state0,hypers,A,xiout,options)
%% PREDICTION Predict new values for GEK
% 
% [xoutn varxoutn predreport] = predict(xi,state0,hypers,A,xiout,options)
%
% Predicts new values at xiout.
%
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% state0 - initial state, size N x (1+ndim)
% hypers - optimized hyperparameters, size 1 x (ndim + 2)
% A - gain matrix, size N(1+ndim) x N(1+ndim)
% xiout - location of predictions, size n x ndim
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% xoutn - predicted normalized mean, size n x 1
% varxoutn - predicted normalized variance, size n x 1
% predreport - << under construction >>
%
% see also GEK, DEFAULTOPTS, NORMALIZE, HYPERESTIMATE, INITIALSTATE,
% PREDICT, DENORMALIZE, GEKPART1, GEKPART2

%% == general ==
tic;
nc = length(xi);
prefactor = nc / (nc-options.nregresscoeff);
[nout ndim] = size(xiout);

%% == check for covariables ==
[n1 m1] = size(xi);
[n2 m2] = size(state0);
if n1 == n2
    covars = 'no';
else
    covars = 'yes';
end

if options.debug > 0
    if strcmp(options.predicttype,'predict')
        fprintf('gek -> predicting new values ...\n')
    else
        fprintf('gek -> predicting new values and gradients ...\n')
    end
end

%% == split hypers ==
% errnf = nan(n2,1);
% if n1 == n2
%     errnf = hypers(1) * errn;
% else
%     errnf(1:n1) = hypers(1) * errn(1:n1);
%     errnf(n1+1:end) = hypers(2) * errn(n1+1:end);
% end

range = hypers(3:end);

%% == prediction ==
% predict values for xiout:
% If necessary, this is looped, because the system memory
% limits the array size of the correlation matrix 'b'.
% The maxArraySize might be set in options; by default
% it is set according to a system memory check in defaultops.m
% if strcmp(options.estvar,'yes')
%     maxrows = floor( options.maxArraySize / nout / nout );
% else
%     maxrows = floor( options.maxArraySize / nout / ndim );
% end
% nloops = ceil(nout/maxrows);
% xoutn = zeros(nout,1); % preallocate
% varxoutn = zeros(nout,1); % preallocate
% for i = 1:nloops
%     from = (i-1)*maxrows+1;
%     to = min(i*maxrows,nout);
%     b = corrmatrix(xiout(from:to,:),xi,hypers,options.predicttype,covars,options); % b := P * H'
% 
%     % conditional mean
%     xoutn(from:to) = b*state0; % state0 := A^-1 * y := (R + H*P*H')^-1 * y
%     
%     % conditional variance
%     if strcmp(options.estvar,'yes')
%         p = corrmatrix(xiout(from:to,:),xiout(from:to,:),hypers,options.predicttype,'no',options);
%         varxoutn(from:to) = prefactor * diag( p - b*(A\b') );
%     end
% end


% conditional mean
b = corrmatrix(xiout,xi,range,options.predicttype,covars,options); % b := P * H'
xoutn = b*state0; % state0 := A^-1 * y := (R + H*P*H')^-1 * y

% conditional variance
varxoutn = 0*xoutn + NaN;
if strcmp(options.estvar,'yes') && nc > 1
    b = corrmatrix(xiout,xi,range,options.predicttype,covars,options); % b := P * H'
    %p = corrmatrix(xiout,xiout,hypers,'predict','no',options);
    if strcmp(options.predicttype,'predict')
        diagp = ones(nout,1);
    elseif strcmp(options.predicttype,'predictplus')
        diagp = ones(nout,1)*[1 hypers.^-2];
        diagp = reshape(diagp,numel(diagp),1);
    end
    varxoutn = prefactor * ( diagp - diag(b*(A\b')) );
elseif strcmp(options.estvar,'eigen')
    b = corrmatrix(xiout,xi,range,'predict',covars,options); % b := P * H'
    [V D] = eig(A);
    %varxoutn0 = 1 - diag( b*(A\b');
    for k = 1:size(D,1)
        Dmin = D(k,k);
        Vmin = V(:,k);
        scale =  Vmin'*Vmin;
        bmin = 0*b;
        for i = 1:size(b,1);
            bmin(i,:) = b(i,:) * Vmin / scale * Vmin';
        end
        energy(:,k) = diag( bmin * (1/Dmin) * bmin');
    end
    varxoutn = 1 - sum(energy,2);
end

%% == report ==
predreport.cpu = toc;
%predreport.memoryloops = nloops;
predreport.memoryloops = NaN;
end