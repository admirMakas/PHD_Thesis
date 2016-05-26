function [ Y_hat,MSE,dY_hat_dx ] = gekPred( model,testpoints )

% Load Values
R = model.R;
F = model.F;
beta = model.beta;
Var = model.Var;
%MLE = model.MLE;
n = model.n;
%nprime = model.nprime;
S = model.S;
Y = model.Y;
%Ymin = model.ymin;
%xi = model.INPUTS.xi;
%y = model.INPUTS.y;

[row, col] = size(testpoints);

testpoints = (testpoints - kron(model.Norm.NormX(1,:),ones(row,1)) ) ./ kron(model.Norm.NormX(2,:),ones(row,1));
 %kron(NormX(1,:),ones(rows,1)

% % prod
% Corrgauss = @(xi,xj) exp(-model.theta.*(xi - xj).^2);
% Corrgauss_xi = @(xi,xj)  -2*model.theta.*(xi-xj).*exp(-model.theta.*(xi-xj).^2);
% Corrgauss_xi_xi = @(xi,xj)  2*model.theta.* exp(-model.theta.*(xi-xj).^2).*(2* model.theta.*(xi-xj).^2-1);


Y_hat = zeros(row,1);
MSE = zeros(row,1);
%EI = zeros(length(testpoints),1);
dY_hat_dx = zeros(row,model.k);

for i = 1:row
%     r = zeros(length(S),1);
%     dr = r;
%     for j = 1:length(S)
%         if j<=n
%             r(j) = Corrgauss( model.Norm.S(j,:),testpoints(i,:));
%             dr(j) = Corrgauss_xi( model.Norm.S(j,:),testpoints(i,:) );
%         else
%             r(j) = Corrgauss_xi( model.Norm.S(j,:),testpoints(i,:) );
%             dr(j) = Corrgauss_xi_xi( model.Norm.S(j,:),testpoints(i,:) );
%         end
%     end

    [r,dr] = findRs(model.theta, model.Norm.S,  testpoints(i,:),model.n, model.k); 
    
    Y_hat(i) = beta + transpose(r)*(R\(model.Norm.Y-beta*F));
    MSE(i) = Var*(1.0 - transpose(r)*(R\r) + (transpose(r)*(R\F) - 1)^2 / (transpose(F)*(R\F)));
    dY_hat_dx(i,:) = -1*transpose(dr)*(R\(model.Norm.Y-beta*F));
    %EI(i) = (Ymin - Y_hat(i)) * (0.5 + 0.5 *erf( (Ymin - Y_hat(i))/ sqrt(Var*2) ) + (sqrt(Var)/(sqrt(2*pi)))*exp(-1*( (Ymin-Y_hat(i))^2 /(2*Var))));
    
end

Y_hat = Y_hat.*model.Norm.NormY(2,:) + model.Norm.NormY(1,:);
dY_hat_dx = dY_hat_dx.*kron(model.Norm.NormX(1,:),ones(row,1))./kron(model.Norm.NormY(1,:),ones(row,col));

% Y_hat = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* Y_hat;
% %MSE = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* MSE;

end

function [r,dr] = findRs(theta, S, testpoints,n, k) 

Corrgauss = @(xi,xj) exp(-sum(theta.*(xi - xj).^2));
Corrgauss_xi = @(xi,xj,l) exp(-sum(theta.*(xi - xj).^2)).*(-2*theta(l)...
    .*(xi(1,l) - xj(1,l)));
%==========================================================================
Corrgauss_xi2 = @(xi,xj,l) exp(-sum(theta.*(xi - xj).^2)).*...
    (2*theta(l) - 4*theta(l).^2*(xi(1,l) - xj(1,l)).^2);
Corrgauss_xi_xj = @(xi,xj,l,m) exp(-sum(theta.*(xi - xj).^2)).*...
    (-4*theta(l).*theta(m).*(xi(1,l) - xj(1,l)) .* (xi(1,m) - xj(1,m)));

%==========================================================================
% Y_hat = zeros(length(testpoints),1);
% MSE = zeros(length(testpoints),1);
% EI = zeros(length(testpoints),1);
%==========================================================================
i = 1;
%for i = 1:length(testpoints(:,1))
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
    dr = reshape(dr,[(k+1)*n,k]);
    % Y_hat(i) = beta + r'*inv(R)*(Y-beta*F);
    % ======================================================================
    %Y_hat(i) = beta + r'*(R\(Y-beta*F));
    %MSE(i) = Var*(1.0 - r'*inv(R)*r + (r'*inv(R)*F - 1)^2 / (F'*inv(R)*F));
    %=======================================================================
    
    %EI(i) = (Ymin - Y_hat(i)) * (0.5 + 0.5 *erf( (Ymin - Y_hat(i))/ sqrt(Var*2) ) + (sqrt(Var)/(sqrt(2*pi)))*exp(-1*( (Ymin-Y_hat(i))^2 /(2*Var))));
    
%end



end