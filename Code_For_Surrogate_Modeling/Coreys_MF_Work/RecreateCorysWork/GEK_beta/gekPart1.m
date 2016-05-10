function gekmodel = ...
    gekPart1(xi,y,erry,grad,errgrad,options)
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
