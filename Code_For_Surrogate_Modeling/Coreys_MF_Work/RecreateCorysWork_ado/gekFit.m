function [ model ] = gekFit( xi,xigrads, y, grad)

[n,dim] = size(xi);
[nprime,dimprime] = size(xigrads);

if dim == 1;
    [xi,index] = sort(xi);
    y = y(index);
    [xigrads,index] = sort(xigrads);
    grad = grad(index);
end

[Norm] = NORM( xi, y, grad);

S = [xi;xigrads];
Y = [y;reshape(grad,[nprime*dimprime,1])];

% Input section
model.INPUTS.xi = xi;
model.INPUTS.xigrads = xigrads;
model.INPUTS.y = y;
model.INPUTS.grad = grad;
model.Norm = Norm;

% Other information
model.n = n;
model.nprime = nprime;
model.S = S;
model.Y = Y;
model.ymin = min(y);
model.k = dim;

% regression F
model.F = [ones(n,1);zeros(nprime*dimprime,1)];
F = model.F;

% Theta bounds for optimization
UTheta=ones(1,dim).*1.3010;%1.3010; was 2
LTheta=ones(1,dim).*-1; % was -3

% Genetic Alg.
[thetas,~]=gaF(@(x) MLE(x, n, nprime, Norm.xi, Norm.Y, F),1,[],[],[],[], LTheta,UTheta);
% was S now xi 

model.theta = 10.^thetas;

[model.MLE,model.R,model.beta,model.Var] = MLE(thetas, n, nprime, Norm.xi, Norm.Y, F);


end

%CREATE FUNCTION TO GET LIKELIHOOD
function [MLE_val,R,beta,Var] = MLE(x, n, nprime, S, Y, F)

thetas=10.^x;

% R  = createCovMatrix(thetas,S,n,nprime);
[ R ] = BuildRdot( S,thetas );

% [L,U]=lu(R);
% beta = ((F'*(U\(L\F)))\F')*(U\(L\Y));
% Var = (1/(n+nprime)) * (Y-beta*F)'*(U\(L\(Y-beta*F)));
% MLE_val = real(-(n+nprime)*log(Var) - log( det(L)*det(U) ));
% MLE_val = -1*MLE_val;
%Old Code==================================================================
%F = [ones(n,1);zeros(nprime,1)];
beta = ((F'*(R\F))\F')*(R\Y);
Var = (1/(n+nprime)) * (Y-beta*F)'*(R\(Y-beta*F));
MLE_val = real(-(n+nprime)*log(Var) - log( det(R) ));
MLE_val = -1*MLE_val;
%==========================================================================
end

function [Norm] = NORM( xi, y, grad)

NormX = [mean(xi);std(xi)];
NormY = [mean(y);std(y)];

%ones(2,2)*NormX(1,:)
[rows, col] = size(xi);
xi = ( xi- kron(NormX(1,:),ones(rows,1) ) ) ./ kron(NormX(2,:),ones(rows,1));



y= ( y-NormY(1,:) ) ./ NormY(2,:);
grad = grad.*kron(NormX(2,:),ones(rows,1))./( NormY(2,:).*ones(rows,col));


[nprime,dimprime] = size(grad);

S = [xi;xi];
Y = [y;reshape(grad,[nprime*dimprime,1])];


Norm = struct('NormX',NormX, 'NormY',NormY, 'xi',xi, ...
    'y',y, 'grad',grad, 'S',S,'Y',Y);

end





% function [ R ] = createCovMatrix(theta,S,n,nprime)
% 
% Corrgauss = @(xi,xj) exp(-theta*(xi - xj)^2);
% Corrgauss_xi = @(xi,xj) exp(-theta*(xi-xj)^2)*(-2*theta*(xi-xj));
% Corrgauss_xj = @(xi,xj) exp(-theta*(xi-xj)^2)*(2*theta*(xi-xj));
% Corrgauss_xi_xj_ieqj = @(xi,xj) exp(-theta*(xi-xj)^2)*(-2*theta*(xi-xj))*(2*theta*(xi-xj)) + ...
%     exp(-theta*(xi-xj)^2)*2*theta;
% %Corrgauss_xi_xj_iNOTeqj = @(xi,xj) -4*theta
% 
% 
% TopLeft = zeros(n,n);
% 
% for i = 1:n
%     for j = 1:n
%         TopLeft(i,j) =  Corrgauss( S(i) , S(j) );
%     end
% end
% 
% if nprime == 0
%     
%     R = TopLeft;
% else
%     
%     TopRight = zeros(n,nprime);
%     for i = 1:n
%         for j = 1:nprime
%             TopRight(i,j) =  Corrgauss_xj( S(i) , S(n+j) );
%         end
%     end
%     
%     bottomLeft = zeros(nprime,n);
%     for i = 1:n
%         for j = 1:nprime
%             bottomLeft(j,i) = Corrgauss_xi( S(n+j) , S(i) );
%         end
%     end
%     
%     bottomRight = zeros(nprime,nprime);
%     for i = 1:nprime
%         for j = 1:nprime
% %             if i == j
%                 bottomRight(j,i) = Corrgauss_xi_xj_ieqj( S(n+j) , S(n+i) );
% %             elseif
% %                 
% %             end
%         end
%     end
%      R = [TopLeft,TopRight; bottomLeft, bottomRight];
%     %R = [TopLeft,TopRight; TopRight, bottomRight];
% end
% 
% 
% end





% % Plots of theta 
% thetas = linspace(-3, 2,1000)';
% likhoods = thetas;
% F = [ones(n,1);zeros(nprime*dimprime,1)];
% 
% for i =1:length(thetas)
%    likhoods(i) =  MLE(thetas(i), n, nprime, xi, Y, F);
% end
% 
% figure
% plot(thetas,likhoods)
