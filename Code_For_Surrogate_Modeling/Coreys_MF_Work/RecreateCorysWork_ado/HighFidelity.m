function [ y, dy_dx ] = HighFidelity( x )
%y = zeros(length(x),1);
%dy_dx = y;
y = ((6.*x-2).^2).*sin(12.*x-4);
dy_dx =  (12.*(6*x-2).^2).*cos(4-12.*x)-12.*(6.*x-2).*sin(4-12.*x);

end

