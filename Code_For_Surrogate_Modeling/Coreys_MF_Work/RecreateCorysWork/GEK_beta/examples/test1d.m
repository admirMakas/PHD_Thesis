clc, clear, close all

addpath ..
addpath ../functions

%% fix random states
rand('state',0)
randn('state',0)

%% samples
N = 6; % number of samples ("cfd runs")
xi = lhsdesign(N,1); % location of samples ("design of experiment")
            % so xi might denote a number of angles of attack for which you
            % are going to compute the drag
xiout = linspace(0,1,1e3)'; % output grid for plotting
            % so xiout are the AoA for which you want to 'predict' the drag
            % with the response surface

%% black box solver
% this can be replaced by the cfd code
func = @(xi) sin(4*pi*xi);
%func = @(xi) xi.^2 + xi;
%func = @(xi) 100 * (xi-0.5).^2 .* exp(-32*(xi-0.5).^2)

% accuracy of the cfd code
erry = 0; % measurement error for value ("discretization error of cfd code")
errgr = 0; % measurement error for gradient (usually larger than erry)

% evaluate black box ("cfd code") at locations xi
% % note that we use a weird notation where x denotes the true value (for
% % example the drag), while y denotes an observation (for example a drag
% % measurement in the wind tunnel or a drag computed from cfd)
y = func(xi) + erry*randn(N,1); % values + noise

% here we use "complex step derivatives" to obtain gradients, in cfd code
% we use "adjoint approach"
h = 1e-9;
grad = imag(func(complex(xi,h)))/h + errgr*randn(N,1); % gradients + noise

%% Kriging / GEK options
% % options you might want to use, see 'defaultopts.m' for more

% debugging info:
options.debug = 2; % gives output for debugging, set to 0 for no output

% maximum likelihood settings:
% options.hyperinit = [1 1 .1];
% options.hyperspace = [0 0 1];
% options.hyperest = 'brute';
% options.brutesize = 1e3;

%% Kriging
% most importantly, the gradient input is left empty ('') to use Kriging
% without gradient information
[xout0 varxout gradout vargradout report] = ...
    gek(xi,y,erry*ones(size(y)),'','',xiout,options);

% the RMS Error of the response surface
rmseKriging = sqrt(mean((xout0-func(xiout)).^2))

%% GEK
% now, we do include gradient information
[xout varxout gradout vargradout report] = ...
    gek(xi,y,erry*ones(size(y)),grad,errgr*ones(size(grad)),xiout,options);

% the RMS Error of the response surface
rmseGEK = sqrt(mean((xout-func(xiout)).^2))


%% figs
d = 0.05;
figure(1)
plot(xiout,func(xiout),'k-',...
    xi,y,'bo',...
    xiout,xout0,'g-',xiout,xout,'r-',...
    xiout,xout+sqrt(varxout),'r:',xiout,xout-sqrt(varxout),'r:',...
    [xi-d xi+d]',[y-d*grad y+d*grad]','b-')
legend('test function','sampled values','Kriging','GEK','GEK error band')
xlabel('\xi'), ylabel('x (and y)'), title('response surface')
axis([0 1 -1.5 1.5])
