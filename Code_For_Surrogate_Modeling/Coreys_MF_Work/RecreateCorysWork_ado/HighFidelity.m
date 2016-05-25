function [ y, dy_dx ] = HighFidelity( x )

[y, dy_dx ] = twoD(x);

end

function [y, dy_dx ] = oneD(x)

y = ((6.*x-2).^2).*sin(12.*x-4);
dy_dx =  (12.*(6*x-2).^2).*cos(4-12.*x)-12.*(6.*x-2).*sin(4-12.*x);

end


function [y, dy_dx ] = twoD( x )

    [ y,dy_dx ] = TwoD_High( x );

end