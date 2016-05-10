function P = corrmatrix(xi1,xi2,hypers,type,covars,options)
%% CORRMATRIX Corrmatrix (for GEK)
%
% P = corrmatrix(xi1,xi2,hypers,type,covars,options)
%
% Correlation matrix for Gradient-Enhanced Kriging
%
% see also HYPERESTIMATE, INITIALSTATE, PREDICT

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
[nxi1 ndim] = size(xi1);
[nxi2 ndim] = size(xi2);
                
%% == lag ==

% check memory and preallocate lag
if options.debug > 0
    mag = nxi1*nxi2*ndim; % required memory
    if mag > options.maxArraySize
        warning('memory might be insufficient to construct lag matrix')
    end
end
H = zeros(nxi1,nxi2,ndim); Hn = H; % preallocate

% slow implementation
% tic
% for i = 1:nxi1
%     for j = 1:nxi2
%         for k = 1:ndim
%             H(i,j,k) = xi1(i,k) - xi2(j,k);
%             Hn(i,j,k) = ( xi1(i,k) - xi2(j,k) ) / hypers(k) ;
%         end
%     end
% end
% toc

% fast implementation
%tic
for k = 1:ndim
    XI1 = xi1(:,k)*ones(1,nxi2);
    XI2 = ones(nxi1,1)*xi2(:,k)';
    H(:,:,k) = XI1 - XI2;
    Hn(:,:,k) = (XI1 - XI2)/hypers(k);
end
%toc

%% == correlation matrix ==

% determine number of rows of matrix blocks
% for the analysis we need the full matrix
% for the prediction we only need the top row
if strcmp(type,'analyze') && strcmp(covars,'yes')
    rows = ndim + 1; % construct full matrix: all rows of blocks
elseif strcmp(type,'predictplus')
    rows = ndim + 1;
else
    rows = 1; % construct only top row of blocks
end

if strcmp(covars,'yes')
    cols = ndim + 1;
elseif strcmp(covars,'no')
    cols = 1;
end

% check memory and preallocate correlation matrix
if options.debug > 0
    mag = nxi1*rows*nxi2*(ndim+1); % required memory
    if mag > options.maxArraySize
        warning('memory might be insufficient to construct correlation matrix')
    end
end
P = zeros(nxi1*rows,nxi2*cols);

% precompute matrices
Hn2 = sum(Hn.^2,3);
EXP = exp(-0.5*Hn2);

% fill blocks
for i = 1:rows
    for j = 1:cols
        if i == 1 && j == 1 % no derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                EXP;
        elseif i == 1 % first derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                H(:,:,j-1) .* hypers(j-1)^-2 .* EXP;
        elseif j == 1 % first derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                - H(:,:,i-1) .* hypers(i-1)^-2 .* EXP;
        elseif i == j % second derivative
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                hypers(i-1)^-2 * EXP - ...
                H(:,:,i-1).^2 .* hypers(i-1)^-4 .* EXP;
        else % two partial derivatives
            P(nxi1*(i-1)+1:nxi1*i,nxi2*(j-1)+1:nxi2*j) = ...
                - H(:,:,i-1) .* hypers(i-1)^-2 .* ...
                H(:,:,j-1) .* hypers(j-1)^-2 .* ...
                EXP;
        end
    end
end


