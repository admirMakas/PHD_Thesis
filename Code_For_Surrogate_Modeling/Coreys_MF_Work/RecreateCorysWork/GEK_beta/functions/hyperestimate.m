function [hypers hypreport] = ...
    hyperestimate(xi,yn,errn,options)
%% HYPERESTIMATE Estimate hyperparameters for GEK
% 
% [hypers hypreport] = ...
%    hyperestimate(xi,yn,errn,options)
%
% Optimizes the hyperparameters for GEK, using a Maximum Likelihood 
% Estimate (MLE). 
%
% Note that in the code there are ndim + 2 hyperparameters, i.e. one
% erryfactor, one errgradfactor, and ndim correlation ranges. By default,
% GEK does not optimize for the err factors, and optimizes for all
% correlation ranges.
%
% Optimizing the hyperparameters is often the most expensive part of GEK in
% terms of CPU time, so it might pay to provide good initiall guesses and
% consider the relevant options, such as:
%
% options.debug = 0; prints optimization results <default>
% options.debug = 1; prints optimization results
% options.debug = 2; prints & plots optimization results
%
% options.hyperest = 'fmin'; uses fminsearch optimization <default>
% options.hyperest = 'brute'; uses brute force optimization
%
% options.hyperinit = [1 1 0.3 0.4]; this is an example setting for the
%  initial guess of the hyperparameters for a 2-d case, where erryfact = 1, 
%  errgradfact = 1, correlation range1 = 0.3, and range2 = 0.4
%  <default setting depends on N and on domain size>
%
% options.hyperspace = [0 0 1 1]; 2-d example, optimizes only for both 
%  correlation ranges <default>
% options.hyperspace = [0 1 0 1]; 2-example, optimizes only for the
%  errgradfact and for correlation range2
%
% options.brutesize = 256; # samples in brute optimization <default>
% 
% N - number of observations
% n - number of predictions
% ndim - number of spatial dimensions
%
% == INPUT ==
% xi - location of observations, size N x ndim
% yn - normalized data vector, size N x (1+ndim)
% errn - normalized error vector, size N x (1+ndim)
% options - GEK options, refer to defaultopts
%
% == OUTPUT ==
% hypers - optimized hyperparameters, size 1 x (ndim + 2)
% hypreport - << under construction >>
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
[ny ndim] = size(xi);
[nyn dum] = size(yn);

if nyn > ny
    gradients = 'yes';
else
    gradients = 'no';
end

tic

%% == estimate hypers ==
switch options.hyperest
    case 'none'
        if options.debug > 0
            fprintf('gek -> fixed hyperparameters \n')
        end
        hypers = options.hyperinit;
        erFact = 'empty';
    case 'brute'
        if options.debug > 0
            fprintf('gek -> estimating hyperparameters (brute force) ... \n')
        end
        hypspace = options.hyperspace;
        hypinit = options.hyperinit;
        hypid = find(hypspace~=0);
        hypdim = length(hypid);
        hyp0 = lhsdesign(options.brutesize,hypdim,'iter',256);
        hyp0 = (options.brutefactor).^(2*hyp0-1);
        hyp = repmat(hypinit(hypid),options.brutesize,1) .* hyp0;
        loghyp = log(hyp);
        for i = 1:options.brutesize
            negloglike(i) = ...
                goalfunction(loghyp(i,:),hypid,hypinit,...
                    xi,yn,errn,options.goalfunc,gradients);
        end
        [m im] = min(negloglike);
        hypers = hypinit;
        hypers(hypid) = exp(loghyp(im,:));
        if hypdim == 1 && options.debug > 1
            figure(10)
            loglog(exp(loghyp),exp(-negloglike),'b.')
            hold on 
            loglog(exp(loghyp(im)),exp(-m),'ko')
            hold off
            xlabel(['hyper ' num2str(hypid)])
            title('likelihood')
        elseif hypdim == 2 && options.debug > 1
            figure(10)
            scatter(log10(exp(loghyp(:,1))),log10(exp(loghyp(:,2))),...
                3,exp(-negloglike))
            hold on
            scatter(log10(exp(loghyp(im,1))),log10(exp(loghyp(im,2))),...
                40,'ko','filled')
            hold off
            xlabel(['log10 hyper ' num2str(hypid(1))])
            ylabel(['log10 hyper ' num2str(hypid(2))])
            colorbar
            title('likelihood')
        end
    case 'fmin'
        hypspace = options.hyperspace;
        hypinit = options.hyperinit;
        hypid = find(hypspace~=0);
        [loghypfmin,fval,exitflag,output] = ...
            fminsearch(@(loghyp)goalfunction(loghyp,hypid,hypinit,...
                    xi,yn,errn,options.goalfunc,gradients),...
                        log(hypinit(hypid)),options.fminopts);
        hypers = hypinit;
        hypers(hypid) = exp(loghypfmin);
    otherwise
        error('this type of minimization of the goalfunction is not supported')
end

if options.debug > 0
    fprintf(['     erryfact    : ' num2str(hypers(1)) '\n'])
    fprintf(['     errgradfact : ' num2str(hypers(2)) '\n'])
    fprintf(['     range (the) : ' num2str(hypers(3:end)) '\n'])
end

%% == report ==
hypreport.cpu = toc;
