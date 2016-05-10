function [yn errn norming normreport] = ...
    normalize(xi,y,erry,grad,errgrad,options)
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
