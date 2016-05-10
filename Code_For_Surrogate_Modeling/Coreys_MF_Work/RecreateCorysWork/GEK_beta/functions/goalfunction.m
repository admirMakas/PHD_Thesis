function measure = goalfunction(loghyper,hypid,hypinit,...
                    xi,yn,errn,type,gradients)
%% GOALFUNCTION Goalfunction (for GEK)
%
% measure = goalfunction(loghyper,hypid,hypinit,xi,yn,errn,type,gradients)
%
% Goal function for optimization of hyperparameters. Calculates negative
% log likelihood.
%
% see also HYPERESTIMATE

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

mleopts.debug = 0;

N = size(xi,1);
N2 = size(yn,1);

%% reconstruct hyper
fullhyper = hypinit;
fullhyper(hypid) = exp(loghyper);

if N == N2
    erryfactor = fullhyper(1);
else
    erryfactor = fullhyper(1);
    errgradfactor = fullhyper(2);
end
range = fullhyper(3:end);

errnf = nan(N2,1);
if N == N2
    errnf = erryfactor * errn;
else
    errnf(1:N) = erryfactor * errn(1:N);
    errnf(N+1:end) = errgradfactor * errn(N+1:end);
end

P = corrmatrix(xi,xi,range,'analyze',gradients,mleopts);             % prior
R = diag(errnf.^2);                       
A = R + P;

switch type
    case 'mle'
        measure = sum(log(eig(A))) + yn'*(A\yn);
    case 'mle_efficient'
        % more efficient, less robust
        [U p] = chol(A);
        if 0 == 0
           % this is more efficient
           measure = sum(log(eig(U))) + yn'*(U'\(U\yn));
        else
           % this is more robust
           measure = sum(log(eig(A))) + yn'*(A\yn);
        end
    case 'mle_selected'
        % ???? this might improve results for small N ????
        f = ones(size(yn,1),1);
        if size(yn,1) > size(xi,1)
            f(size(xi,1)+1:end) = 0; % do we want t
        end
        measure = sum(log(eig(A))) + (f.*yn)'*(A\yn);
    otherwise
        warning('this type of goal function is not supported')
end
