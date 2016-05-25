function [ y, dy_dx ] = LowFidelity( x )
% The cheap model is a modified version of the high fidelity model for this
% equation
 [y, dy_dx ] = twoD( x );


end

function [y, dy_dx ] = oneD(x)

[ fe, dfe_dx ] = HighFidelity( x );

A = 0.5; B = 10; C = -5;
y =  A.*fe + (B.*(x-0.5)) + C; 
dy_dx =  A*dfe_dx + B;

end


function [y, dy_dx ] = twoD( x )

    [ y,dy_dx ] = TwoD_Low( x );

end