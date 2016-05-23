function [ model ] = gekFit( xi,xigrads, y, grad)

[xi,index] = sort(xi);
y = y(index);
[xigrads,index] = sort(xigrads);
grad = grad(index);

S = [xi;xigrads];
Y = [y;grad];

[n, k] = size(xi);
nprime = k*n;

model.INPUTS.xi = xi;
model.INPUTS.xigrads = xigrads;
model.INPUTS.y = y;
model.INPUTS.grad = grad;

model.F = [ones(n,1);zeros(nprime,1)];
    F=model.F;

UTheta=ones(1,1).*1.5;
LTheta=ones(1,1).*-1;

[thetas,NegLnLikelihood]=...
ga(@(x) MLE(x, n, k, nprime, S, Y, F),k,[],[],[],[], LTheta,UTheta);

model.theta = 10.^thetas;

model.R  = createCovMatrix(model.theta,S,n,k);
    R=model.R;
    
%Go with LU decomposition for now
[L,U]=lu(R);
model.beta = ((F'*(U\(L\F)))\F')*(U\(L\Y));
    beta=model.beta;
model.Var = (1/(n+nprime)) * (Y-beta*F)'*(U\(L\(Y-beta*F)));
    Var=model.Var;
model.MLE = (n+nprime)*log(Var) + log( det(L)*det(U) );

model.n = n;
model.k = k;
model.nprime = nprime;
model.S = S;
model.Y = Y;
model.ymin = min(y);

end

%CREATE FUNCTION TO GET LIKELIHOOD
function [MLE_val] = MLE(x, n, k, nprime, xi, Y, F)

thetas=10.^x;

R  = createCovMatrix(thetas,xi,n,k);

[L,U]=lu(R);
beta = ((F'*(U\(L\F)))\F')*(U\(L\Y));
Var = (1/(n+nprime)) * (Y-beta*F)'*(U\(L\(Y-beta*F)));
MLE_val = (real(-(n+nprime)*log(Var) - log( det(L)*det(U) )))*-1;

%Old Code==================================================================
%F = [ones(n,1);zeros(nprime,1)];
% beta = ((F'*(R\F))\F')*(R\Y);
% Var = (1/(n+nprime)) * (Y-beta*F)'*(R\(Y-beta*F));
% MLE_val = real(-(n+nprime)*log(Var) - log( det(R) ));
%==========================================================================
end

% CREATE CORREATION MATRIX
function [ R ] = createCovMatrix(theta,xi,n,k)

% Pre–allocate memory
R1 = zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
    for j=i+1:n
        R1(i,j)=exp(-sum(theta.*(xi(i,:)-xi(j,:)).^2));
    end
end
% Add upper and lower halves and diagonal of ones plus
% small number to reduce ill conditioning
R1=R1+eye(n)+eye(n).*eps;

% Pre–allocate memory
RDot=zeros((k+1)*n, (k+1)*n);
RDot(1:n, 1:n)=R1(1:end, 1:end);
for l=1:k
    % Get first order derivatives
    RDoti = zeros(n, n);
    for i=1:n
        for j=i+1:n
            RDoti(i, j) = 2*theta(l)*(xi(i,l) - xi(j, l))*R1(i,j);
        end
    end
    RDoti=RDoti-RDoti'+zeros(n);
    RDot(1:n, (n*l)+1:(n*l)+n) = RDoti(1:end, 1:end);
    
    % Get second order derivatives
    for m=1:k
    RDoti = zeros(n, n);
        if l==m
            for i=1:n
                for j=i+1:n
                    RDoti(i, j) = (2*theta(l) - 4*theta(l)^2*...
                        (xi(i,l)-xi(j,l))^2)*R1(i,j);
                end
            end
            RDoti=RDoti+(2*theta(l)*eye(n));
            RDot((n*l)+1:(n*l)+n, (n*l)+1:(n*l)+n) = RDoti(1:end, 1:end);
        elseif m>l
            for i=1:n
                for j=i+1:n
                    RDoti(i, j) = -4*theta(l)*theta(m)*(xi(i,l) - xi(j,l))*...
                        (xi(i,m) - xi(j,m))*R1(i,j);
                end
            end
            RDoti=RDoti+RDoti'+zeros(n);
            RDot((n*l)+1:(n*l)+n, (n*m)+1:(n*m)+n) = RDoti(1:end, 1:end);
        end
    end
end
%Assemble RDot matrix
RDot=RDot+(triu(RDot,1))';
R=RDot;
end