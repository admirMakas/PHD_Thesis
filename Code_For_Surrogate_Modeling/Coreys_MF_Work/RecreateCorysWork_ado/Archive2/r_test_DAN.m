clear all
clc

X=[1; 2];
S=[X;X];

[n, k] = size(X);

testpoints = [3]

[row, col] = size(testpoints);

theta=[0.5];

% prod
Corrgauss = @(xi,xj) exp(-theta.*(xi - xj).^2);
Corrgauss_xi = @(xi,xj)  -2*theta.*(xi-xj).*exp(-theta.*(xi-xj).^2);
Corrgauss_xi_xi = @(xi,xj)  2*theta.* exp(-theta.*(xi-xj).^2).*(2*theta.*(xi-xj).^2-1);

for i = 1:row
    r = zeros(length(S),1);
    dr = r;
    for j = 1:length(S)
        if j<=n
            r(j) = Corrgauss( S(j,:),testpoints(i,:));
            dr(j) = Corrgauss_xi( S(j,:),testpoints(i,:) );
        else
            r(j) = Corrgauss_xi( S(j,:),testpoints(i,:) );
            dr(j) = Corrgauss_xi_xi( S(j,:),testpoints(i,:) );
        end
    end
end