function [ indexs ] = ptsCaptured2( center, data, TRS)
Bounds = [center-TRS; center+TRS];
[points, dims] = size(data);
Captured = ones(points,dims);
for  i = 1:points
    for j = 1:dims
        if data(i,j) >= Bounds(1,j) && data(i,j) <= Bounds(2,j)
            Captured(i,j) = 0; 
        end
    end
end
[indexs] = find( sum(Captured,2) == 0 );
end