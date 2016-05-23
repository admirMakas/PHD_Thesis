X=[1 2; 3 4];

% Determine data dimensions
[n, k] = size(X);

theta=[0.5 0.5];

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