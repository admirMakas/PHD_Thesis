function [ Hist ] = TRS_calc( itr, Constraints, Hist )
tol = 0.01;

if itr == 0
    TRS = (Constraints.xMaxs - Constraints.xMins)/4;
    Hist.rho = [Hist.rho; 0];
else
    %       Current          Next                Current          App Next
    rho = ( Hist.yFh(itr) - Hist.yFh(itr+1) ) / ( Hist.yFh(itr) - Hist.obj(itr));
    if itr > 10000 %4
        distances = DISTANCES( Hist.centers(end,:), Hist.centers);
        avgDist = distances(end-2:end-1);
        avgDist = mean(avgDist);
    else
        avgDist = tol+1;
    end
    
    if avgDist < tol
        display('hit')
        TRS = .25*Hist.TRS(itr,:);
    else

        if rho > 0.8
            TRS = 2*Hist.TRS(itr,:); % was 2*TRS This creates errors  .7
        elseif rho > 0.5
            TRS = Hist.TRS(end,:);
        else
            TRS = 0.5*Hist.TRS(itr,:);
        end
        
    end
    
    Hist.rho = [Hist.rho; rho];
end

Hist.TRS = [Hist.TRS; TRS];

end
% Add something like if distance to previous point keeps decreasing or
% average distance, then trust size decreases significantly. THis is
% because as more information is added to the localitites close together
% the surrogate becomes inaccurate.


function distances = DISTANCES( center, data)

center_array = ones(size(data));
for i = 1:length(center)
    center_array(:,i) =  center_array(:,i)*center(i);
end

if length(center) == 1
    
    distances = sqrt((center_array - data).^2 );
    
else
    
    distances = sqrt(sum(transpose( (center_array - data).^2 ))');
    
end

end









% function [ Hist ] = TRS_calcOLD( itr, Constraints, Hist )
% 
% if itr == 0
%     TRS = (Constraints.xMaxs - Constraints.xMins)/2;
%     Hist.rho = [Hist.rho; 0];
% else
%     %       Current          Next                Current          App Next
%     rho = ( Hist.yFh(itr) - Hist.yFh(itr+1) ) / ( Hist.yFh(itr) - Hist.obj(itr));
% 
%         
%         
%     if rho > 0.8
%         TRS = .7*Hist.TRS(itr); % was 2*TRS This creates errors  .7
%     elseif rho > 0.5
%         TRS = Hist.TRS(end);
%     else
%         TRS = 0.5*Hist.TRS(itr);
%     end
% 
%     
%     Hist.rho = [Hist.rho; rho];
% end
% 
% Hist.TRS = [Hist.TRS; TRS];
% 
% end

