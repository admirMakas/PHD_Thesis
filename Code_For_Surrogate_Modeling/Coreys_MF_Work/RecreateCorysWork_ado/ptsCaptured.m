function [ indexs ] = ptsCaptured( center, data, distance)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

center_array = ones(size(data));
for i = 1:length(center)
    center_array(:,i) =  center_array(:,i)*center(i);
end

if length(center) == 1
    
    distances = sqrt((center_array - data).^2 );
    
else
    
    distances = sqrt(sum(transpose( (center_array - data).^2 ))');
    
end

indexs = find(distances < distance);

end

%% Old code 
% rows = length(data);
% indexs = zeros(rows,1);
% indexs = [];
% k = 1;
% for i = 1:rows
%     if ( abs(center - data(i)) < distance ) == 1
%         indexs(i,1) = i;
%     else
%         indexs(k) = i;
%         k = k+1;
%     end
%     % distances(i,1) = ( abs(center - xdata(i)) < distance );
% end
% indexs(indexs) = [];