function [ f,df ] = TwoD_High( x )

[row,col] = size(x);
f = zeros(row,1);
df = zeros(row,col);

for i =1:row
    [f(i),df(i,:)] = rosenbrock(x(i,:));
end


end

function [f,df] = rosenbrock(x)
    %f = sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1:end-1)-1).^2);
    
    
    f = 100*(x(:,2)-x(:,1).^2).^2+(1-x(:,1)).^2;


    
    if nargout > 1
        df = [-400*(x(:,2)-x(:,1)^2)*x(:,1)-2*(1-x(:,1)) ; 200*(x(:,2)-x(:,1)^2)];
%         n = length(x);
%         df(1,2:n-1) = -400*reshape(x(2:n-1),n-2,1).*(reshape(x(3:n),n-2,1)...
%             -reshape(x(2:n-1),n-2,1).^2)+2*(reshape(x(2:n-1),n-2,1)-1)+...
%             200*(reshape(x(2:n-1),n-2,1)-reshape(x(1:n-2),n-2,1).^2);
%         df(1,1) = -400*(x(2)-x(1)^2)*x(1)+2*(x(1)-1);
%         df(1,n) = 200*(x(n)-x(n-1)^2);
    end
end