function [ f,df ] = TwoD_Low( x )

[row,col] = size(x);
f = zeros(row,1);
df = zeros(row,col);

for i =1:row
    [f(i),df(i,:)] = sphere(x(i,:));
end


end

function [f,df] = sphere(x)
    f = sum(100*x.^2);
    if nargout > 1
        df = 200*x;
    end
end