function [ Y_hat,MSE ] = gekPred( model,testpoints )

% Load Values
R = model.R;
F = model.F;
beta = model.beta;
Var = model.Var;
MLE = model.MLE;
n = model.n;
k = model.k;
nprime = model.nprime;
S = model.S;
Y = model.Y;
Ymin = model.ymin;
xi = model.INPUTS.xi;
y = model.INPUTS.y;

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
% 

Corrgauss = @(xi,xj) exp(-sum(model.theta.*(xi - xj).^2));
Corrgauss_xi = @(xi,xj,l) exp(-sum(model.theta.*(xi-xj).^2))*(-2*model.theta(l)*(xi(1,l)-xj(1,l)));

Y_hat = zeros(length(testpoints),1);
MSE = zeros(length(testpoints),1);
EI = zeros(length(testpoints),1);


for i = 1:length(testpoints(:,1))
    %r = zeros(3*n,1);
    r = [];
    
    for j = 1:2*n %1:length(S)
        if j<=n
            r(j,1) = Corrgauss( S(j,:),testpoints(i,:) )';
        elseif j>n
            for l=1:k
                r1(l,1) = Corrgauss_xi( S(j-n,:),testpoints(i,:),l);
            end
            r=[r;r1];
        end
    end
    
    % Y_hat(i) = beta + r'*inv(R)*(Y-beta*F);
    Y_hat(i) = beta + r'*(R\(Y-beta*F));
    MSE(i) = Var*(1.0 - r'*inv(R)*r + (r'*inv(R)*F - 1)^2 / (F'*inv(R)*F));
    
    %EI(i) = (Ymin - Y_hat(i)) * (0.5 + 0.5 *erf( (Ymin - Y_hat(i))/ sqrt(Var*2) ) + (sqrt(Var)/(sqrt(2*pi)))*exp(-1*( (Ymin-Y_hat(i))^2 /(2*Var))));

end

% Y_hat = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* Y_hat;
% %MSE = repmat(model.par.Ysc(1,:),mx,1) + repmat(model.par.Ysc(2,:),mx,1) .* MSE;

end
