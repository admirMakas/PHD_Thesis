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

testpoints = (testpoints - model.Norm.NormX(1,:) ) ./ model.Norm.NormX(2,:);


% prod
Corrgauss = @(xi,xj) exp(-model.theta.*(xi - xj).^2);
Corrgauss_xi = @(xi,xj)  -2*model.theta.*(xi-xj).*exp(-model.theta.*(xi-xj).^2);
Corrgauss_xi_xi = @(xi,xj)  2*model.theta.* exp(-model.theta.*(xi-xj).^2).*(2* model.theta.*(xi-xj).^2-1);


Y_hat = zeros(length(testpoints),1);
MSE = zeros(length(testpoints),1);
EI = zeros(length(testpoints),1);
dY_hat_dx = Y_hat;

for i = 1:length(testpoints)
    r = zeros(length(S),1);
    dr = r;
    for j = 1:length(S)
        if j<=n
            r(j) = Corrgauss( model.Norm.S(j,:),testpoints(i,:));
            dr(j) = Corrgauss_xi( model.Norm.S(j,:),testpoints(i,:) );
        else
            r(j) = Corrgauss_xi( model.Norm.S(j,:),testpoints(i,:) );
            dr(j) = Corrgauss_xi_xi( model.Norm.S(j,:),testpoints(i,:) );
        end
    end
    % Y_hat(i) = beta + r'*inv(R)*(Y-beta*F);
    Y_hat(i) = beta + transpose(r)*(R\(model.Norm.Y-beta*F));
    MSE(i) = Var*(1.0 - transpose(r)*(R\r) + (transpose(r)*(R\F) - 1)^2 / (transpose(F)*(R\F)));
    dY_hat_dx(i) = -1*transpose(dr)*(R\(model.Norm.Y-beta*F));
    %EI(i) = (Ymin - Y_hat(i)) * (0.5 + 0.5 *erf( (Ymin - Y_hat(i))/ sqrt(Var*2) ) + (sqrt(Var)/(sqrt(2*pi)))*exp(-1*( (Ymin-Y_hat(i))^2 /(2*Var))));
    
end

Y_hat = Y_hat.*model.Norm.NormY(2,:) + model.Norm.NormY(1,:);
dY_hat_dx = dY_hat_dx.*model.Norm.NormX(1,:)./model.Norm.NormY(1,:);

% Y_hat = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* Y_hat;
% %MSE = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* MSE;



end
