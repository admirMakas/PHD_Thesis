function options = ...
    defaultopts(options,xi)
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
