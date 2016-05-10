function [xout varxout gradout vargradout denormreport] = ...
    denormalize(xiout,xoutn,varxoutn,norming,options)
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
    vargradout = reshape(vargradout,nx,ndim);
else
    gradout = NaN;
    vargradout = NaN;
end

%% == report ==
denormreport.cpu = toc;
