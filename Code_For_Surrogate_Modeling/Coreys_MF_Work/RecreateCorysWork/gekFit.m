function [ model ] = gekFit( xi,xigrads, y, grad)
tic
n = length(xi);
nprime = length(xigrads);

S = [xi;xigrads];
Y = [y;grad];

model.INPUTS.xi = xi;
model.INPUTS.xigrads = xigrads;
model.INPUTS.y = y;
model.INPUTS.grad = grad;

thetas = linspace(1,10,6000)';
thetas = 10.^thetas;

MLEs = zeros(length(thetas),1);

% USES THE PREDEFINED VECTOR OF THETA VALUES ABOVE RANGING FROM 1-10
% AND CALCULATES ESTIMATES FOR THE MEAN (beta) AND VARIANCE (var)
for Q = 1:length(thetas)
    theta = thetas(Q);
    R  = createCovMatrix(theta,S,n,nprime);
    F = [ones(n,1);zeros(nprime,1)];
    beta = ((F'*(R\F))\F')*(R\Y);
    Var = (1/(n+nprime)) * (Y-beta*F)'*(R\(Y-beta*F));
    MLE = real(-(n+nprime)*log(Var) - log( det(R) ));
    if MLE ==   Inf
        MLEs(Q) = Inf;
    else
        MLEs(Q) = -1*MLE;
    end
end

%figure(10)
%plot(thetas,MLEs)
% BASED ON THE ABOVE LOOP CALCULATIONS BELOW CODE EXTACTS THE THETA
% THAT GIVES MIN MLE ESTIMATE
model.theta = thetas(find( MLEs == min(MLEs)));
if length(model.theta)>1
    model.theta = model.theta(1);
end


model.R  = createCovMatrix(model.theta,S,n,nprime);
model.F = [ones(n,1);zeros(nprime,1)];
model.beta = ((F'*(R\F))\F')*(R\Y);
model.Var = (1/(n+nprime)) * (Y-beta*F)'*(R\(Y-beta*F));
model.MLE = (n+nprime)*log(Var) + log( det(R) );
model.n = n;
model.nprime = nprime;
model.S = S;
model.Y = Y;
model.ymin = min(y);
model.thetas = log10(thetas);
model.MLEs = MLEs;

toc
end

function [ R ] = createCovMatrix(theta,S,n,nprime)

Corrgauss = @(xi,xj) exp(-theta*(xi - xj)^2);
Corrgauss_xi = @(xi,xj) exp(-theta*(xi-xj)^2)*(-2*theta*(xi-xj));
Corrgauss_xj = @(xi,xj) exp(-theta*(xi-xj)^2)*(2*theta*(xi-xj));
Corrgauss_xi_xj = @(xi,xj) exp(-theta*(xi-xj)^2)*(-2*theta*(xi-xj))*(2*theta*(xi-xj)) + ...
    exp(-theta*(xi-xj)^2)*2*theta;

TopLeft = zeros(n,n);

for i = 1:n
    for j = 1:n
        TopLeft(i,j) =  Corrgauss( S(i) , S(j) );
    end
end

if nprime == 0
    
    R = TopLeft;
else
    
    TopRight = zeros(n,nprime);
    for i = 1:n
        for j = 1:nprime
            TopRight(i,j) =  Corrgauss_xj( S(i) , S(n+j) );
        end
    end
    
    bottomLeft = zeros(nprime,n);
    for i = 1:n
        for j = 1:nprime
            bottomLeft(j,i) = Corrgauss_xi( S(n+j) , S(i) );
        end
    end
    
    bottomRight = zeros(nprime,nprime);
    for i = 1:nprime
        for j = 1:nprime
            bottomRight(j,i) = Corrgauss_xi_xj( S(n+j) , S(n+i) );
        end
    end
    R = [TopLeft,TopRight; bottomLeft, bottomRight];
end


end