function [ yhat,dyhat_dx,intercept ] = LinearModel_dan( estPt, point, y, dy_dxi, intercept)

if nargin == 5
    
    yhat = intercept + estPt*dy_dxi';
    dyhat_dx = dy_dxi;
else
    intercept = findIntercept(point, y, dy_dxi);
    yhat = 0;
    dyhat_dx = 0';
end


end

function [intercept] = findIntercept(point, y, gradients)
intercept = y-(point*gradients');

end