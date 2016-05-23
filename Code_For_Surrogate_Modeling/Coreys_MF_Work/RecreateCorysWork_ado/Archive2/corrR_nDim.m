clear all
clc
tic
X=[-0.7071; 0.7071];

% Determine data dimensions
[n, k] = size(X);

theta=[5.5372];

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

%==========================================================================

% Pre–allocate memory
RDot = zeros((k+1)*n, (k+1)*n);
% Build upper half of PsiDot
for l=1:k
    for m=l:k
        if l==1
        % Build upper half of dPsidX
            for i=1:n
                for j=i+1:n
                    RDot(i, m*n+j) = 2*theta(m)*(X(i,m) - X(j, m))*R(i,j);
                end
            end
        % Add upper and lower halves
        RDot(1:n, m*n+1:(m+1)*n) = RDot(1:n, m*n+1:(m+1)*n)...
        - RDot(1:n, m*n+1:(m+1)*n)';
        end 
        
        if m==l
        % Build upper half of d2PsidX ?2
            for i=1:n
                for j=i+1:n
                    RDot(l*n+i, m*n+j)=...
                    (2*theta(l) - 4*theta(l)^2*(X(i,l)-X(j,l))^2)*R(i,j);
                end
            end
        % Add half diagonal
        RDot(l*n+1:(l+1)*n, m*n+1:(m+1)*n)=...
        RDot(l*n+1:(l+1)*n, m*n+1:(m+1)*n) + eye(n).*theta(l);
    
        else
        % Build upper half of d2PsidXdX
            for i=1:n
                for j=i+1:n
                    RDot(l*n+i, m*n+j)=...
                    -4*theta(l)*theta(m)*(X(i,l) - X(j,l))*(X(i,m) - X(j,m))...
                    *R(i,j);
                end
            end
        % Add upper and lower halves
        RDot(l*n+1:(l+1)*n, m*n+1:(m+1)*n)=...
        RDot(l*n+1:(l+1)*n, m*n+1:(m+1)*n) + RDot(l*n+1:(l+1)...
        *n,m*n+1:(m+1)*n)';
        end
    end
end
% Add upper and lower halves to Psi
RDot=[R zeros(n,k*n); zeros(k*n,(k+1)*n)] + RDot + RDot'
toc