function [state0 A initreport] = ...
    initialstate(xi,yn,errn,hypers,options)
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

% 2012, JHS de Baar, TU Delft
%
% Kluyverweg 1 (blue high rise) - Room 10.18
% 2629 HS Delft
% The Netherlands
%
% Phone: +31 15 2782596
% Cell Phone: +31 6 42703414
% e-Mail:  j.h.s.debaar[insert_curly]tudelft.nl
%
% Feel free to contact me with any questions on the application of GEK or
% requests for new or taylormade versions of this code.
%
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>

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
