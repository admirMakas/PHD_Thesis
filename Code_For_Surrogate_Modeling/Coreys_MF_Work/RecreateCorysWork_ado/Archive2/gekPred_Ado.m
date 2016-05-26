%function [ Y_hat,MSE ] = gekPred( model,testpoints )
clear all
clc
tic

X=[-0.16298 0.519863;
    -0.73433 -1.0113;
    1.695957 -0.03165;
    -1.50543 1.552218;
    1.049683 -1.59302];

S=[X;X];

testpoints=[0.620689655172414 1.31034482758621];

theta=[0.523558 0.1];

[n, k] = size(X);

%==========================================================================
% Load Values
% R = model.R;
% F = model.F;
% beta = model.beta;
% Var = model.Var;
% MLE = model.MLE;
% n = model.n;
% nprime = model.nprime;
% S = model.S;
% Y = model.Y;
% Ymin = model.ymin;
% xi = model.INPUTS.xi;
% y = model.INPUTS.y;
%==========================================================================

% % Normalize testpoints
% [m n] = size(xi);  % number of design sites and number of dimensions
% sx = size(testpoints);            % number of trial sites and their dimension
% if  min(sx) == 1 & n > 1 % Single trial point
%     nx = max(sx);
%     if  nx == n
%         mx = 1;  testpoints = testpoints(:).';
%     end
% else
%     mx = sx(1);  nx = sx(2);
% end
% if  nx ~= n
%     error(sprintf('Dimension of trial sites should be %d',n))
% end
% x_original=testpoints;
% testpoints = (testpoints - repmat(model.par.Ssc(1,:),mx,1)) ./ repmat(model.par.Ssc(2,:),mx,1);

%====================================================================================================
% Corrgauss = @(xi,xj) exp(-sum(model.theta.*(xi - xj).^2));
% Corrgauss_xi = @(xi,xj,l,j) exp(-sum(model.theta.*(xi-xj).^2))*(-2*model.theta(l)*(xi(j,l)-xj(j,l)));
%====================================================================================================

Corrgauss = @(xi,xj) exp(-sum(theta.*(xi - xj).^2));
Corrgauss_xi = @(xi,xj,l) exp(-sum(theta.*(xi - xj).^2)).*(-2*theta(l)...
    .*(xi(1,l) - xj(1,l)));
%==========================================================================
% Corrgauss_xi2 = @(xi,xj,l) exp(-sum(theta.*(xi - xj).^2)).*((2*theta(l))...
%     *(-2*theta(l)*(xi(1,l) - xj(1,l)).^2 + 1));
Corrgauss_xi2 = @(xi,xj,l) exp(-sum(theta.*(xi - xj).^2)).*...
    (2*theta(l) - 4*theta(l).^2*(xi(1,l) - xj(1,l)).^2);
Corrgauss_xi_xj = @(xi,xj,l,m) exp(-sum(theta.*(xi - xj).^2)).*...
    (-4*theta(l).*theta(m).*(xi(1,l) - xj(1,l)) .* (xi(1,m) - xj(1,m)));

%==========================================================================
% Y_hat = zeros(length(testpoints),1);
% MSE = zeros(length(testpoints),1);
% EI = zeros(length(testpoints),1);
%==========================================================================

for i = 1:length(testpoints(:,1))
    %r = zeros(3*n,1);
    r = [];
    dr = [];
    %======================================================================
    
    for j = 1:n %1:length(S)
        r(j,1) = Corrgauss( S(j,:),testpoints(i,:) )';
    end
    
    for l=1:k
        for j=1:n
            r1(j,1) = Corrgauss_xi( S(j,:),testpoints(i,:),l);
        end
        r=[r;r1];
    end
    %======================================================================
    for l=1:k
    dr1 = [];    
    RDoti = zeros(n, 1);
        for j=1:n
            RDoti(j,1) = Corrgauss_xi( S(j,:),testpoints(i,:),l);
        end
        dr1=[dr1;RDoti];
        
        for m=1:k
            RDoti = zeros(n, 1);
            if l==m
                for j=1:n
                    RDoti(j,1) = Corrgauss_xi2( S(j,:),testpoints(i,:),l);
                end
                dr1=[dr1;RDoti];
            else
                for j=1:n
                    RDoti(j,1) = Corrgauss_xi_xj( S(j,:),testpoints(i,:),l,m);
                end
                dr1=[dr1;RDoti];
            end
        end
    dr=[dr;dr1];
    end
    r
    dr = reshape(dr,[(k+1)*n,k])
    % Y_hat(i) = beta + r'*inv(R)*(Y-beta*F);
    % ======================================================================
    %Y_hat(i) = beta + r'*(R\(Y-beta*F));
    %MSE(i) = Var*(1.0 - r'*inv(R)*r + (r'*inv(R)*F - 1)^2 / (F'*inv(R)*F));
    %=======================================================================
    
    %EI(i) = (Ymin - Y_hat(i)) * (0.5 + 0.5 *erf( (Ymin - Y_hat(i))/ sqrt(Var*2) ) + (sqrt(Var)/(sqrt(2*pi)))*exp(-1*( (Ymin-Y_hat(i))^2 /(2*Var))));
    
end

% Y_hat = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* Y_hat;
% %MSE = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* MSE;

%end
toc