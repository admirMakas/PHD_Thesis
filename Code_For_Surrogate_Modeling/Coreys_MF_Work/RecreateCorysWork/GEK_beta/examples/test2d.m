clear, close all

addpath ..
addpath ../functions

%% fix random states
rand('state',0)
randn('state',0)

%% samples
N = 20; % number of samples ("cfd runs")
xi = lhsdesign(N,2); % location of samples ("design of experiment")
            % in multi-dimensional, each column of xi represents a
            % dimension
xiout = lhsdesign(2e3,2); % output grid for plotting
            % you could also use a meshgrid, then you can have contour
            % plots, however, kriging can handle any grid

% xiout = 2*lhsdesign(8e3,2)-0.5; % output grid for plotting
            % uncomment the above line, to see what happens with the 
            % response surface when you move outside the sampling region

%% black box solver
% this can be replaced by the cfd code
func = @(xi) sin(4*pi*sum(xi,2))+cos(3*pi*xi(:,1));

% accuracy of the cfd code
erry = .1; % measurement error for value ("discretization error of cfd code")
errgr = .25; % measurement error for gradient (usually larger than erry)

% evaluate black box ("cfd code") at locations xi
y = func(xi) + erry*randn(N,1); % values + noise

% here we use "complex step derivatives" to obtain gradients, in cfd code
% we use "adjoint approach"
ndim = size(xi,2);
grad = nan(N,ndim);
h = 1e-9;
for i = 1:ndim
    xih = xi;
    xih(:,i) = complex(xih(:,i),h);
    grad(:,i) = imag(func(xih))/h + errgr*randn(N,1); % gradients + noise
end

%% Kriging / GEK options
% options you might want to use, see 'defaultopts.m' for more

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

%% plotting
exact = func(xiout);
cmin = min(exact);
cmax = max(exact);

figure(1)
scatter(xiout(:,1),xiout(:,2),3,exact), colorbar
title('true function'), xlabel('\xi_1'), ylabel('\xi_2')

figure(2)
scatter(xiout(:,1),xiout(:,2),3,xout0), colorbar
hold on, scatter(xi(:,1),xi(:,2),50,'ko','filled'), hold off
title('Kriging response surface'), xlabel('\xi_1'), ylabel('\xi_2')

figure(3)
scatter(xiout(:,1),xiout(:,2),3,xout), colorbar
hold on, scatter(xi(:,1),xi(:,2),50,'ko','filled'), hold off
title('GEK response surface'), xlabel('\xi_1'), ylabel('\xi_2')