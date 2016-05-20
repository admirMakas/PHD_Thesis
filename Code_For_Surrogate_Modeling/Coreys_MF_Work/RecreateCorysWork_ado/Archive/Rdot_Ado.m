clear all
clc
tic
X=[1 2 3; 3 4 5];

% Determine data dimensions
[n, k] = size(X);

theta=[0.5 0.5 0.5];

% Pre–allocate memory
R=zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
    for j=i+1:n
        R(i,j)=exp(-sum(theta.*(X(i,:)-X(j,:)).^2));
    end
end
% Add upper and lower halves and diagonal of ones plus
% small number to reduce ill conditioning
R=R+R'+eye(n)+eye(n).*eps;

% Pre–allocate memory
RDotR = [];, RDotC = [];
RDot=zeros((k+1)*n, (k+1)*n);
% Build upper half of PsiDot
for l=1:k
    RDoti = zeros(n, n);
    for i=1:n
        for j=i+1:n
            RDoti(i, j) = -2*theta(l)*(X(i,l) - X(j, l))*R(i,j);
        end
    end
    RDoti=RDoti-RDoti'+zeros(n);
    RDotC=[RDotC;RDoti];
    RDotR=cat(2, R, RDotC');
end

RDot2=[];
for l=1:k
    RDot2i=[];
    for m=1:k
    RDoti = zeros(n, n);
        if l==m
            for i=1:n
                for j=i+1:n
                    RDoti(i, j) = (2*theta(l) - 4*theta(l)^2*...
                        (X(i,l)-X(j,l))^2)*R(i,j);
                end
            end
            RDoti=RDoti+RDoti'+eye(n);
            RDot2i=[RDot2i;RDoti];
        else
            for i=1:n
                for j=i+1:n
                    RDoti(i, j) = -4*theta(l)*theta(m)*(X(i,l) - X(j,l))*...
                        (X(i,m) - X(j,m))*R(i,j);
                end
            end
            RDoti=RDoti+RDoti'+zeros(n);
            RDot2i=[RDot2i;RDoti];
        end
    end
    RDot2=cat(1, RDot2, RDot2i');
end

RDot(1:n, 1:(k+1)*n) = RDotR(1:end, 1:end);
RDot(n+1:end, 1:n) = RDotC(1:end, 1:end);
RDot(n+1:end, n+1:end) = RDot2(1:end, 1:end);
toc