function [ Hist ] = nDimensionalMF( TYPE, center, Bounds, options, PLOT,gridSize)
ind = 1; ConvergenceCheck = 1000;

[ Constraints ] = problemConstraints( Bounds(1,:), Bounds(2,:) );
[ Hist ] = MultiFidelityHist( );
[ Hist ] = TRS_calc( 0, Constraints, Hist );
Hist.centers = [Hist.centers; center];
[ Hist ] = EvaluateFidelities( Hist.centers(1,:), Hist );

if strcmp(PLOT,'YES')
    if length(center) == 1 % 1D
        testpoints = linspace(Constraints.xMins,Constraints.xMaxs,gridSize)';
    elseif length(center) == 2 % 2D
        if gridSize > 200
            gridSize = 30;
        end
        testpoints = gridsamp(Bounds, gridSize);
    else
        display('Good luck plotting this.')
        PLOT = 'NO'
        
    end
    
end

while ConvergenceCheck > 0.0000001
    [ Hist ] = AdditiveAndMultiplicative( ind, Hist );
    [ pointsCaped ] = ptsCaptured2( Hist.centers(ind,:), Hist.centers,Hist.TRS(ind,:)./2 );
    
    if length(pointsCaped) < 2
        display('linear')
        lin = 1;
        Hist.w1 = [Hist.w1; 0.5]; Hist.w2 = [Hist.w2; 0.5];
        
        [ ~,~ ,interceptAdd ] = LinearModel_dan( [], Hist.centers(ind,:), ...
            Hist.yAdd(pointsCaped), Hist.dAdd_dx(pointsCaped,:));
        
        [ ~,~,interceptMult ] = LinearModel_dan( [], Hist.centers(ind,:), ...
            Hist.yMult(pointsCaped), Hist.dMult_dx(pointsCaped,:));
        
        intercepts = [interceptAdd; interceptMult];   
        
    else
        display('GEK')
        lin = 0;
        % Addative Model
        xi = Hist.centers(pointsCaped,:); xigrads = xi;
        y = Hist.yAdd(pointsCaped); grad = Hist.dAdd_dx(pointsCaped,:);
        [ modelAdd ] = gekFit( xi,xigrads, y, grad);
        % [ modelAdd ] = gekFit( xi,[], y, []);
        
        % Multiplicitive Model
        y = Hist.yMult(pointsCaped); grad = Hist.dMult_dx(pointsCaped,:);
        [ modelMult ] = gekFit( xi,xigrads, y, grad);
        
        y = Hist.yFh(pointsCaped); grad = Hist.dFh_dx(pointsCaped,:);
        [ modelHIGH ] = gekFit( xi,xigrads, y, grad);
        
        % [ modelMult ] = gekFit( xi,[], y, []);
        Hist.HighModels = [Hist.HighModels;modelHIGH];
        Hist.AddModels = [Hist.AddModels;modelAdd];
        Hist.MultModels = [Hist.MultModels;modelMult];
        
        % Prediction of each model
        ADD = @(x) gekPred( modelAdd,x );
        MULT = @(x) gekPred( modelMult,x );
        
        % ERROR function
        Error = @(y,yhat) sqrt( sum( (y - yhat).^2)/length(y) );
        
        % Finding the weights of the models
        PHI = @(N,sigma) (1/(2*pi*sigma^2))^(N/2)*exp(-N/2);
        
        x_removeCaped = Hist.centers; x_removeCaped(pointsCaped,:) = [];
        Fh_removeCaped = Hist.yFh; Fh_removeCaped(pointsCaped) = [];
        Fl_removeCaped = Hist.yFl; Fl_removeCaped(pointsCaped) = [];
        
        if isempty(x_removeCaped)
            Hist.w1 = [Hist.w1; 0.5];
            Hist.w2 = [Hist.w2; 0.5];
        else
            
            if strcmp(TYPE,'DAN')
                sigmaAdd = Error(Fh_removeCaped, Fl_removeCaped+ADD(x_removeCaped));
                sigmaMult = Error(Fh_removeCaped, Fl_removeCaped.*MULT(x_removeCaped));
                
                phiAdd = PHI(length(Fh_removeCaped),sigmaAdd);
                phiMult = PHI(length(Fh_removeCaped),sigmaMult);
                
                wi = ((Hist.w1(ind-1)*phiAdd)/((Hist.w1(ind-1)*phiAdd)+(Hist.w2(ind-1)*phiMult)));
                Hist.w1 = [Hist.w1; wi];
                wi = ((Hist.w2(ind-1)*phiAdd)/((Hist.w1(ind-1)*phiAdd)+(Hist.w2(ind-1)*phiMult)));
                Hist.w2 = [Hist.w2; wi];
            else
                
                [ pointsCapedwR ] = ptsCaptured2( Hist.centers(ind,:), x_removeCaped, 1.5*Hist.TRS(ind,:)./2 );
                if isempty(pointsCapedwR)
                    display('Cory 50/50')
                    Hist.w1 = [Hist.w1; 0.5];
                    Hist.w2 = [Hist.w2; 0.5];
                else
                    display('Cory Weights')
                    x_removeCaped = x_removeCaped(pointsCapedwR,:);
                    Fh_removeCaped = Fh_removeCaped(pointsCapedwR);
                    Fl_removeCaped = Fl_removeCaped(pointsCapedwR);
                    
                    sigmaAdd = Error(Fh_removeCaped, Fl_removeCaped+ADD(x_removeCaped));
                    sigmaMult = Error(Fh_removeCaped, Fl_removeCaped.*MULT(x_removeCaped));
                    
                    phiAdd = PHI(length(Fh_removeCaped),sigmaAdd);
                    phiMult = PHI(length(Fh_removeCaped),sigmaMult);
                    
                    wi = ((Hist.w1(ind-1)*phiAdd)/((Hist.w1(ind-1)*phiAdd)+(Hist.w2(ind-1)*phiMult)));
                    Hist.w1 = [Hist.w1; wi];
                    wi = ((Hist.w2(ind-1)*phiAdd)/((Hist.w1(ind-1)*phiAdd)+(Hist.w2(ind-1)*phiMult)));
                    Hist.w2 = [Hist.w2; wi];
                end
                
            end
        end
        
    end
    
    [ LB, UB, Hist ] = AdjustBounds( Hist, Constraints );    

    if lin == 1
        % [F,grad] = HighFidelityApprox_linear( [-1,-1], Hist, intercepts, pointsCaped, Hist.w1(end), Hist.w2(end) )
        [center,objectiveValue,~,OUTPUTCON] = fmincon(...
            @(x)HighFidelityApprox_linear(x,Hist,intercepts,pointsCaped,Hist.w1(end), ...
            Hist.w2(end)),Hist.centers(end,:),[],[],[],[],LB, UB,[],options);
    else
        %[F,grad] = HighFidelityApprox([-.5,-.5],modelAdd,modelMult,Hist.w1(ind), Hist.w2(ind)),Hist.centers(end,:);
        
        [center,objectiveValue,~,OUTPUTCON] = fmincon(...
            @(x)HighFidelityApprox(x,modelAdd,modelMult,Hist.w1(ind), ...
            Hist.w2(ind)),Hist.centers(end,:),[],[],[],[],LB, UB,[],options);
    end
    Hist.optzLcount = Hist.optzLcount + OUTPUTCON.funcCount;
    Hist.obj = [Hist.obj; objectiveValue];
    Hist.centers = [Hist.centers; center];
    
    [ Hist ] = EvaluateFidelities( Hist.centers(end,:), Hist );
    
    if strcmp(PLOT,'YES')
        if length(Constraints.xMins) == 1
            
            if lin == 1
                
                oneDplot(testpoints, Hist, LB, UB, intercepts, pointsCaped)
            else
                
                oneDplot(testpoints, Hist, LB, UB, [], [], modelAdd, modelMult)
            end
            
            
        else
            if lin == 1
                
                twoDplot(testpoints, Hist, LB, UB, intercepts, pointsCaped)
            else
                
                twoDplot(testpoints, Hist, LB, UB, [], [], modelAdd, modelMult)
            end
            
            
            
        end
        
    end
    
    [ Hist ] = TRS_calc( ind, Constraints, Hist ); 
    ConvergenceCheck = (abs(Hist.yFh(ind) - Hist.yFh(ind+1))) / abs(Hist.yFh(ind));
    
    % ConvergenceCheck = (abs(Hist.yFh(ind) - Hist.yFh(ind+1))) / abs(Hist.yFh(ind)-objectiveValue);
    ind = ind + 1;

end

Hist.Results = [length(Hist.yFh), length(Hist.yFh)+Hist.optzLcount, Hist.centers(end)];

end

function [] = oneDplot(testpoints, Hist, LB, UB, intercepts, pointsCaped, modelAdd, modelMult)

figure, hold on
plot(testpoints,LowFidelity(testpoints))
plot(testpoints,HighFidelity(testpoints))
legend('Y_{Cheap}','Y_{Expensive}')
plot([LB,UB],[-10,-10],'g','linewidth',2)
scatter(Hist.centers(end-1),-10,50,'fill','k')
if nargin == 6
    plot(testpoints,HighFidelityApprox_linear(testpoints,Hist,...
        intercepts,pointsCaped,Hist.w1(end),Hist.w2(end)),...
        'k','LineWidth',2)
else
    plot(testpoints,HighFidelityApprox(testpoints,modelAdd,...
        modelMult,Hist.w1(end), Hist.w2(end)),'k','LineWidth',2)
end

scatter(Hist.centers, Hist.yFh,200,'o')
scatter(Hist.centers(end,:), Hist.yFh(end),50,'fill','k')
scatter(Hist.centers, Hist.yFl,200,'o')
scatter(Hist.centers(end,:), Hist.yFl(end),50,'fill','k')
axis([0  1.1 -10 20])

end



function [] = twoDplot(testpoints, Hist, LB, UB, intercepts, pointsCaped, modelAdd, modelMult)
Lr = LowFidelity(testpoints);
Hr = HighFidelity(testpoints);
gridSize = 30;

X1 = reshape(testpoints(:,1),gridSize,gridSize);
X2 = reshape(testpoints(:,2),gridSize,gridSize);
ResponseL = reshape(Lr, size(X1));
ResponseH = reshape(Hr, size(X1));

% Plotting the window size
figure(2), hold on
rectangle('Position',[LB UB(1)-LB(1) UB(2)-LB(2) ])
text(UB(1)+abs(UB(1)*.1),LB(2),int2str(length(Hist.centers(1:end-1,1))))
scatter(Hist.centers(end-1,1),Hist.centers(end-1,2),50,'fill','k')
scatter(Hist.centers(end,1),Hist.centers(end,2),50,'fill','r')
scatter(1,1,50,'fill','b')
axis([-2,2,-2,2])

% Plotting the approxmation to check stability
% figure, hold on
% mesh(X1,X2,ResponseA)
% view([-20,25])

% History of evaluations
figure
subplot(2,2,1), hold on
contour(X1,X2,ResponseL,30)
scatter(Hist.centers(1:end-1,1),Hist.centers(1:end-1,2),50,'fill','k')
scatter(Hist.centers(end,1),Hist.centers(end,2),50,'fill','r')
scatter(1,1,50,'fill','b')
title('Low Fidelity')

subplot(2,2,2), hold on
contour(X1,X2,ResponseH,50)
scatter(Hist.centers(1:end-1,1),Hist.centers(1:end-1,2),50,'fill','k')
scatter(Hist.centers(end,1),Hist.centers(end,2),50,'fill','r')
scatter(1,1,50,'fill','b')
title('High Fidelity')

if nargin == 6
    Ar = HighFidelityApprox_linear(testpoints,Hist,intercepts,pointsCaped,Hist.w1(end),Hist.w2(end));

else
    Ar = HighFidelityApprox(testpoints,modelAdd,modelMult,Hist.w1(end), Hist.w2(end));
end
ResponseA = reshape(Ar, size(X1));

if nargin == 6
    subplot(2,2,[3,4]), hold on
    contour(X1,X2,ResponseA,30)
    rectangle('Position',[LB UB(1)-LB(1) UB(2)-LB(2) ])
    scatter(Hist.centers(1:end-1,1),Hist.centers(1:end-1,2),50,'fill','k')
    scatter(1,1,50,'fill','b')
    title('Approximation')
else
    
    HHr = gekPred(Hist.HighModels(end),testpoints);
    ResponseHHr = reshape(HHr, size(X1));
    
    
    subplot(2,2,3), hold on
    contour(X1,X2,ResponseHHr,50)
    rectangle('Position',[LB UB(1)-LB(1) UB(2)-LB(2) ])
    scatter(Hist.centers(1:end-1,1),Hist.centers(1:end-1,2),50,'fill','k')
    scatter(1,1,50,'fill','b')
    title('High Fidelity Approximation')
    
    subplot(2,2,4), hold on
    contour(X1,X2,ResponseA,50)
    rectangle('Position',[LB UB(1)-LB(1) UB(2)-LB(2) ])
    scatter(Hist.centers(1:end-1,1),Hist.centers(1:end-1,2),50,'fill','k')
    scatter(1,1,50,'fill','b')
    title('Cory''s Method Approximation')
end


end





