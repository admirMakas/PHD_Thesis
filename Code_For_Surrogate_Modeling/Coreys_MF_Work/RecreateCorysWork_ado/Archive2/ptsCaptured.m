function [ distances ] = ptsCaptured( center, xdata, distance)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
rows = length(xdata);
distances = zeros(rows,1);
indexs = [];
k = 1;
for i = 1:rows
    if ( abs(center - xdata(i)) < distance ) == 1
        distances(i,1) = i;
    else
        indexs(k) = i;
        k = k+1;
    end
    % distances(i,1) = ( abs(center - xdata(i)) < distance );
end
distances(indexs) = [];

end