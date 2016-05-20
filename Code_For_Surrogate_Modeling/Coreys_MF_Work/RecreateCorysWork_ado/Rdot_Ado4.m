clear all
clc
tic
X=[1 2 3 4 5 6 7; 3 4 5 6 7 8 9; 6 7 8 9 10 11 12];

% Determine data dimensions
[n, k] = size(X);

theta=[1 0.75 0.5 0.4 0.3 0.25 0.2];

% Preľallocate memory
R = zeros(n,n);
% Build upper half of correlation matrix
for i=1:n
    for j=i+1:n
        R(i,j)=exp(-sum(theta.*(X(i,:)-X(j,:)).^2));
    end
end
% Add upper and lower halves and diagonal of ones plus
% small number to reduce ill conditioning
R=R+eye(n)+eye(n).*eps;

% Preľallocate memory
RDot=zeros((k+1)*n, (k+1)*n);
RDot(1:n, 1:n)=R(1:end, 1:end);
for l=1:k
    % Get first order derivatives
    RDoti = zeros(n, n);
    for i=1:n
        for j=i+1:n
            RDoti(i, j) = 2*theta(l)*(X(i,l) - X(j, l))*R(i,j);
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
                        (X(i,l)-X(j,l))^2)*R(i,j);
                end
            end
            RDoti=RDoti+(2*theta(l)*eye(n));
            RDot((n*l)+1:(n*l)+n, (n*l)+1:(n*l)+n) = RDoti(1:end, 1:end);
        elseif m>l
            for i=1:n
                for j=i+1:n
                    RDoti(i, j) = -4*theta(l)*theta(m)*(X(i,l) - X(j,l))*...
                        (X(i,m) - X(j,m))*R(i,j);
                end
            end
            RDoti=RDoti+RDoti'+zeros(n);
            RDot((n*l)+1:(n*l)+n, (n*m)+1:(n*m)+n) = RDoti(1:end, 1:end);
        end
    end
end
%Assemble RDot matrix
RDot=RDot+(triu(RDot,1))'
toc