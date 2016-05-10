function [xout varxout gradout vargradout report] = ...
    gek(xi,y,erry,grad,errgrad,xiout,options)
%% GEK Gradient-Enhanced Kriging
% 
% [xout varxout gradout vargradout report] = ...
%    gek(xi,y,erry,grad,errgrad,xiout,options)
%
% We consider a Qauntity of Interest (QoI). For N observations y of the QoI 
% at the locations xi, Kriging predicts the mean xout and variation varxout
% of the QoI at the output locations xiout. (For Kriging, supply '' for
% grad and errgrad.)
%
% For the same QoI, if we observe N values y and N x ndim gradients grad,
% GEK predicts xout and varxout, as well as the gradients gradout and
% gradient variations vargradout at the output locations xiout.
%
% In many cases, the data is available before xiout has been determined. In
% these cases, one might want to use gekPart1 and gekPart2.
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
% xout - predicted mean of QoI, size n x 1
% varxout - predicted variance of QoI, size n x 1
% gradout - predicted gradients, size n x ndim
% vargradout - predicted variance of gradients, size n x ndim
% report - << under construction >>
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
    fprintf('\n=== gek ===\n')
end

%% === normalize ===
[yn errn norming normreport] = ...
    normalize(xi,y,erry,grad,errgrad,options);

%% === estimate hyperparameters ===
% this is still incomplete
[hypers hypreport] = ...
    hyperestimate(xi,yn,errn,options);

%% === find initial state ===
[state0 A initreport] = ...
    initialstate(xi,yn,errn,hypers,options);

%% === make prediction ===
[xoutn varxoutn predreport] = ...
    prediction(xi,state0,hypers,A,xiout,options);

%% === denormalize prediction ===
[xout varxout gradout vargradout denormreport] = ...
    denormalize(xiout,xoutn,varxoutn,norming,options);

%% === create final report ===
report = 'to do';
%    createreport(normreport,hypreport,initreport,predreport,denormreport);
