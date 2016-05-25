function [ LB, UB, Hist ] = AdjustBounds( Hist, Constraints )


Bounds = [Hist.centers(end,:)-Hist.TRS(end,:)./2; Hist.centers(end,:)+Hist.TRS(end,:)./2];

RealBounds = [Constraints.xMins; Constraints.xMaxs];
[points, dims] = size(RealBounds);

for j = 1:dims
    
    % if its less than the real lower bound
    if Bounds(1,j) < RealBounds(1,j)
            Bounds(1,j) = RealBounds(1,j);
    end
    
    % if its greater than the real upper bound
    if Bounds(2,j) > RealBounds(2,j)
            Bounds(2,j) = RealBounds(2,j);
    end
    
end
LB = Bounds(1,:);
UB = Bounds(2,:);
%Hist.TRS(end,:) = UB - LB;

end

