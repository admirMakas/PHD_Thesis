function [xoutn varxoutn predreport] = ...
    prediction(xi,state0,hypers,A,xiout,options)
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
