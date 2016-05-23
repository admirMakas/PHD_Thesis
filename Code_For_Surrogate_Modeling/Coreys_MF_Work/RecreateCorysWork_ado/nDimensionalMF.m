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
        testpoints = gridsamp(Bounds, gridSize);
    else
        display('Good luck plotting this.')
        return
    end
    
end


while ConvergenceCheck > 0.01
    [ Hist ] = AdditiveAndMultiplicative( ind, Hist );
    distance = Hist.TRS(ind,:)/2;
    [ pointsCaped ] = ptsCaptured( Hist.centers(ind,:), Hist.centers, distance);
    if length(pointsCaped) < 2
        display('linear')
        lin = 1;
        Hist.w1 = [Hist.w1; 0.5]; Hist.w2 = [Hist.w2; 0.5];
        
        [ ~,~ ,interceptAdd ] = LinearModel_dan( [], Hist.centers(ind,:), ...
            Hist.yAdd(pointsCaped), Hist.dAdd_dx(pointsCaped));
        
        [ ~,~,interceptMult ] = LinearModel_dan( [], Hist.centers(ind,:), ...
            Hist.yMult(pointsCaped), Hist.dMult_dx(pointsCaped));
        
        intercepts = [interceptAdd; interceptMult];   
        
    else
        display('GEK')
        lin = 0;
        % Addative Model
        xi = Hist.centers(pointsCaped,:); xigrads = xi;
        y = Hist.yAdd(pointsCaped); grad = Hist.dAdd_dx(pointsCaped);
        [ modelAdd ] = gekFit( xi,xigrads, y, grad);
        % [ modelAdd ] = gekFit( xi,[], y, []);
        
        % Multiplicitive Model
        y = Hist.yMult(pointsCaped); grad = Hist.dMult_dx(pointsCaped);
        [ modelMult ] = gekFit( xi,xigrads, y, grad);
        % [ modelMult ] = gekFit( xi,[], y, []);
        
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
                [ pointsCapedwR ] = ptsCaptured( Hist.centers(end,:), x_removeCaped, (1.5*Hist.TRS(ind))/2);
                
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
    
    LB = Hist.centers(ind,:)-distance;
    if LB < Constraints.xMins
        LB = Constraints.xMins;
        Hist.TRS(end) = UB-LB;
    end
    UB = Hist.centers(ind,:)+distance;
    if UB > Constraints.xMaxs
        UB = Constraints.xMaxs;
        Hist.TRS(end) = UB-LB;
    end
    %Hist.TRS(end) = UB-LB;
    
    if lin == 1
        [center,objectiveValue,~,OUTPUTCON] = fmincon(...
            @(x)HighFidelityApprox_linear(x,Hist,intercepts,pointsCaped,Hist.w1(end), ...
            Hist.w2(end)),Hist.centers(end,:),[],[],[],[],LB, UB,[],options);
    else
        [center,objectiveValue,~,OUTPUTCON] = fmincon(...
            @(x)HighFidelityApprox(x,modelAdd,modelMult,Hist.w1(ind), ...
            Hist.w2(ind)),Hist.centers(end,:),[],[],[],[],LB, UB,[],options);
    end
    Hist.optzLcount = Hist.optzLcount + OUTPUTCON.funcCount;
    Hist.obj = [Hist.obj; objectiveValue];
    Hist.centers = [Hist.centers; center];
    
    [ Hist ] = EvaluateFidelities( Hist.centers(end,:), Hist );
    
    if strcmp(PLOT,'YES')
        % Check if 1d or 2d too
        figure, hold on
        plot(testpoints,LowFidelity(testpoints))
        plot(testpoints,HighFidelity(testpoints))
        legend('Y_{Cheap}','Y_{Expensive}')
        plot([LB,UB],[-10,-10],'g','linewidth',2)
        scatter(Hist.centers(ind),-10,50,'fill','k')
        if lin == 1
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
    
    [ Hist ] = TRS_calc( ind, Constraints, Hist ); 
    ConvergenceCheck = (abs(Hist.yFh(ind) - Hist.yFh(ind+1))) / abs(Hist.yFh(ind));
    ind = ind + 1;

end

Hist.Results = [length(Hist.yFh), length(Hist.yFh)+Hist.optzLcount, Hist.centers(end)];

end


