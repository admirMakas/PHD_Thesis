function [xout varxout gradout vargradout report] = ...
    gekPart2(gekmodel,xiout)
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
    denormalize(gekmodel.N,xiout,xoutn,...
        varxoutn,gekmodel.norming,gekmodel.options);

%% === create final report ===
report = 'to do';
%    createreport(normreport,hypreport,initreport,predreport,denormreport);
