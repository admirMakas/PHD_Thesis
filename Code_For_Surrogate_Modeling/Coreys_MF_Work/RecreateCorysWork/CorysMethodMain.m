%clear all
close all
clc

%% 1-D
fe = @(x) ((6.*x-2).^2).*sin(12.*x-4); % expensive
dfe_dx = @(x) (12.*(6*x-2).^2).*cos(4-12.*x)-12.*(6.*x-2).*sin(4-12.*x);
A = 0.5; B = 10; C = -5;
fc = @(x) A.*fe(x) + (B.*(x-0.5)) + C; % cheap 
dfc_dx = @(x) A*dfe_dx(x) + B;

%% 1-D Samples
xMin = 0; xMax = 1.1;
testpoints = linspace(xMin,xMax,1000);
%Xc = linspace(0,1,11)'; Rc = fc(Xc);
%Xe = [0, 0.4, 0.6, 1]'; Re = fe(Xe); % These points most coincide with Xc
%Xc = Xe; Rc = fc(Xc);

%% 
% Initialize hold arrays and weights
x_array = []; yAdd_array = []; yMult_array = []; dAdd_array = []; dMult_array = [];
yFh_array = []; yFl_array = []; dFh_array = []; dFl_array = [];
TRS_array = []; rho_array = []; obj_hist = [];
w1 = 0; w2 = 0;
% Step 1
TRS = (xMax - xMin)/2;
TRS_array = [TRS_array; TRS];
center = .6;
x_array = [x_array; center];

% Step 2
Fh = fe(center); Fl = fc(center);
dFh_dx = dfe_dx(center); dFl_dx = dfc_dx(center);

% holds
yFh_array = [yFh_array; Fh]; yFl_array = [yFl_array; Fl]; 
dFh_array = [dFh_array; dFh_dx]; dFl_array = [dFl_array; dFl_dx];
ConvergenceCheck = 1000;
ii = 1;
while ConvergenceCheck > 0.01

%for kjkj = 1:7
%% Plot for each run
figure, hold on
% scatter(Xc, Rc), scatter(Xe, Re)
plot(testpoints,fc(testpoints))
plot(testpoints,fe(testpoints))
legend('Y_{Cheap}','Y_{Expensive}') % 'Cheap Samples','Expensive Samples', 

% Step 3
Add = Fh - Fl;
dAdd_dx = dFh_dx - dFl_dx;
Mult = Fh./Fl;
dMult_dx = (Fl.*dFl_dx - Fh.*dFh_dx)./Fl.^2;

% holds
yAdd_array = [yAdd_array; Add]; yMult_array = [yMult_array; Mult]; 
dAdd_array = [dAdd_array; dAdd_dx]; dMult_array = [dMult_array; dMult_dx];


% Hidden step: check number of samples within model range: if too low
% construct a first order approximation model, else gek models

xdata = x_array; distance = TRS/2;
[ pointsCaped ] = ptsCaptured( center, xdata, distance);
if length(pointsCaped) < 2
    display('linear')
    ADD = @(x) (dAdd_array(pointsCaped).*(x-center)) + yAdd_array(pointsCaped);
    MULT = @(x) (dMult_array(pointsCaped).*(x-center)) + yMult_array(pointsCaped);
    w1 = 0.5;
    w2 = 1-w1;
    fhat_h = @(x) w1.*(fc(x) + ADD(x) ) + w2.*( fc(x).*MULT(x) );
    
else
    display('GEK')
    [ pointsCapedwR ] = ptsCaptured( center, xdata, (1.5*TRS)/2);
    
    % Addative Model
    xi = x_array(pointsCaped); xigrads = xi;
    y = yAdd_array(pointsCaped); grad = dAdd_array(pointsCaped);
    [ modelAdd ] = gekFit( xi,xigrads, y, grad);
    % [ modelAdd ] = gekFit( xi,[], y, []);

    % Multiplicitive Model
    y = yMult_array(pointsCaped); grad = dMult_array(pointsCaped);
    [ modelMult ] = gekFit( xi,xigrads, y, grad);
    %[ modelMult ] = gekFit( xi,[], y, []);
    
    % Prediction of each model
    ADD = @(x) gekPred( modelAdd,x );
    MULT = @(x) gekPred( modelMult,x );
    
    % ERROR function
    Error = @(y,yhat) sqrt( sum( (y - yhat).^2)/length(y) ); 
    
    % Finding the weights of the models
    PHI = @(N,sigma) (1/(2*pi*sigma^2))^(N/2)*exp(-N/2);
    
    x_removeCaped = x_array; x_removeCaped(pointsCaped) = [];
    
    if isempty(x_removeCaped)
        sigmaAdd = 0;
        sigmaMult = 0;
    
%         phiAdd = PHI(length(Fh_removeCaped),sigmaAdd);
%         phiMult = PHI(length(Fh_removeCaped),sigmaMult);
    
        w1 = .5;
        w2 = 1-w1;
    else
        Fh_removeCaped = yFh_array; Fh_removeCaped(pointsCaped) = [];
        Fl_removeCaped = yFl_array; Fl_removeCaped(pointsCaped) = [];
        sigmaAdd = Error(Fh_removeCaped, Fl_removeCaped+ADD(x_removeCaped));
        sigmaMult = Error(Fh_removeCaped, Fl_removeCaped.*MULT(x_removeCaped));
    
        phiAdd = PHI(length(Fh_removeCaped),sigmaAdd);
        phiMult = PHI(length(Fh_removeCaped),sigmaMult);
    
        w1 = (w1*phiAdd)/( (w1*phiAdd) + (w2*phiMult) );
        w2 = (w2*phiMult)/( (w1*phiAdd) + (w2*phiMult) );
        % w2 = 1-w1;
    end
    
    
    fhat_h = @(x) w1.*(fc(x)' + ADD(x) ) + w2.*( fc(x)'.*MULT(x) );

end

plot(testpoints,fhat_h(testpoints),'k','LineWidth',2)
scatter(x_array, yFh_array,200,'o')
scatter(x_array, yFl_array,200,'o')


LB = [center-TRS/2]; 
if LB < xMin
    LB = xMin;
end
UB = [center+TRS/2]; 
if UB > xMax
    UB = xMax;
end
plot([LB,UB],[-10,-10],'g','linewidth',2)

A = []; B = []; Aeq = []; Beq = [];

options = optimset('Algorithm','interior-point','GradObj','off',...
    'GradConstr','off','TolX',1e-6,'Display','iter', 'PlotFcns',...
    @optimplotfval );

[center,objectiveValue] = fmincon(@(x)fhat_h(x),center,A,B,Aeq,Beq,LB,...
UB,[],options);
obj_hist = [obj_hist; objectiveValue];
x_array = [x_array; center]; % New evaluation site
Fh = fe(center); % new expensive model evaluation
Fl = fc(center); % new cheap model evaluation, should be able to extract this from the optz
dFh_dx = dfe_dx(center); % grad values
dFl_dx = dfc_dx(center); % grad values

yFh_array = [yFh_array; Fh]; yFl_array = [yFl_array; Fl]; 
dFh_array = [dFh_array; dFh_dx]; dFl_array = [dFl_array; dFl_dx];

% maybe a differnt way to do this calculation...
rho = ( yFh_array(ii) - yFh_array(ii+1) ) / ( yFh_array(ii) - objectiveValue );

if rho > 0.8
    TRS = 2*TRS; % was 2*TRS This creates errors  .7
elseif rho > 0.5
    
else
    TRS = 0.5*TRS;
end
TRS_array = [TRS_array;TRS];
rho_array = [rho_array; rho];
ConvergenceCheck = (abs(yFh_array(ii) - yFh_array(ii+1))) / abs(yFh_array(ii));
ii = ii + 1;

end
