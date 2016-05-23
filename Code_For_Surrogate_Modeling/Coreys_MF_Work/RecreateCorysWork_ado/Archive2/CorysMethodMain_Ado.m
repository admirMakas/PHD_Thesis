%clear all
close all
clc

%% 1-D, define high and low fidelity functions
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

% Step 1===================================================================
% DEFINE INITIAL TRS AS HALF OF THE ENTIRE DESIGN SPACE
% NEXT LINE APPENDS REVISED TRS VALUES AS ITERATIONS CONTINUE
TRS = (xMax - xMin)/2;
TRS_array = [TRS_array; TRS];

% CENTER IS APPROXIATELY AT THE MID POINT
% NEXT LINE APPENDS REVISED CENTER VALUES AS ITERATIONS CONTINUE
center = 0.6;
x_array = [x_array; center];

% Step 2===================================================================
% GET INITIAL VALUES OF CHEAP AND EXPENSIVE FUNCTION USING CENTER POINT
% ALSO GET DERIVATIVE INFO
Fh = fe(center); Fl = fc(center);
dFh_dx = dfe_dx(center); dFl_dx = dfc_dx(center);

% holds
% APPEND INITIAL VALUES FOR CHEAP AND EXPENSIVE FUNCTIONS INTO THE
% DESIGNATE STORAGE ARRAY
yFh_array = [yFh_array; Fh]; yFl_array = [yFl_array; Fl]; 
dFh_array = [dFh_array; dFh_dx]; dFl_array = [dFl_array; dFl_dx];

%USED TO DETERMINE IF ALGORITHM CONVERGED TO A SOLUTION
ConvergenceCheck = 1000;

% WHILE LOOP FOR REMAINING STEPS IS NEEDED.
ii = 1;
while ConvergenceCheck > 0.001

%for kjkj = 1:7
%% PLOT CHEAP AND EXPENSIVE FUNCTIONS FOR EACH RUN
figure, hold on
% scatter(Xc, Rc), scatter(Xe, Re)
plot(testpoints,fc(testpoints))
plot(testpoints,fe(testpoints))
legend('Y_{Cheap}','Y_{Expensive}') % 'Cheap Samples','Expensive Samples', 

% Step 3===================================================================
% DEFINE ADDITIVE AND MULTIPLICATIVE TERMS AND DERIVATIVES
Add = Fh - Fl;
dAdd_dx = dFh_dx - dFl_dx;
Mult = Fh./Fl;
dMult_dx = (Fl.*dFh_dx - Fh.*dFl_dx)./Fl.^2;
%dMult_dx = dFh_dx./Fl - (dFl_dx.*Fh)/Fl.^2;

% holds
% APPEND VALUES FOR ADDITIVE AND MULTIPLICATIVE TERMS INTO STORAGE
% MATRICES
yAdd_array = [yAdd_array; Add]; yMult_array = [yMult_array; Mult]; 
dAdd_array = [dAdd_array; dAdd_dx]; dMult_array = [dMult_array; dMult_dx];


% Hidden step: check number of samples within model range: if too low
% construct a first order approximation model, else gek models

xdata = x_array; distance = TRS/2;

% FUNCTION 'ptsCaptured' IDENTIFIES WHICH POINTS ARE TO BE KEPT
% FOR ANALYSIS BASED ON THE TRUST REGION SIZE
[ pointsCaped ] = ptsCaptured( center, xdata, distance);

if length(pointsCaped) < 2
    display('linear')
    ADD = @(x) (dAdd_array(pointsCaped).*(x-center)) + yAdd_array(pointsCaped);
    MULT = @(x) (dMult_array(pointsCaped).*(x-center)) + yMult_array(pointsCaped);
    w1 = 0.5;
    w2 = 1-w1;
    
    % HYBRID FUNCTION
    fhat_h = @(x) w1.*(fc(x) + ADD(x) ) + w2.*( fc(x).*MULT(x) );
    
else
    display('GEK')
    % FUNCTION 'ptsCaptured' IDENTIFIES WHICH POINTS ARE TO BE KEPT
    % FOR ANALYSIS BASED ON THE TRUST REGION SIZE
    [ pointsCapedwR ] = ptsCaptured( center, xdata, (1.5*TRS)/2);
    
    % Additive Model=======================================================
    % VECTORS 'xi' AND 'xigrads' USED TO DEFINE POINTS USED FOR KRIGING
    xi = x_array(pointsCaped); xigrads = xi;
    
    % VECTORS 'y' AND 'grad' USED TO DEFINE OUTPUT POINTS USED FOR GRADIENT
    % ENHANCED KRIGING
    y = yAdd_array(pointsCaped); grad = dAdd_array(pointsCaped);
    
    % FUNCTION CALL TO CREATE ADDITIVE SURROGATE MODEL
    [ modelAdd ] = gekFit(xi, xigrads, y, grad);
    %[ modelAdd ] = gekFit(xi, xigrads, yFh_array, dFh_array);
    % [ modelAdd ] = gekFit( xi,[], y, []);

    % Multiplicitive Model=================================================
    
    % VECTORS 'y' AND 'grad' USED TO DEFINE OUTPUT POINTS USED FOR GRADIENT
    % ENHANCED KRIGING
    y = yMult_array(pointsCaped); grad = dMult_array(pointsCaped);
    
    % FUNCTION CALL TO CREATE MULTIPLICATIVE SURROGATE MODEL
    [ modelMult ] = gekFit(xi, xigrads, y, grad);
    %[ modelMult ] = gekFit(xi, xigrads, yFh_array, dFh_array);
    %[ modelMult ] = gekFit( xi,[], y, []);
    
    % Prediction of each model
    ADD = @(x) gekPred( modelAdd,x );
    MULT = @(x) gekPred( modelMult,x );
    
    % ERROR function
    Error = @(y,yhat) sqrt( sum( (y - yhat).^2)/length(y) ); 
    
    % Finding the weights of the models
    PHI = @(N,sigma) (1/(2*pi*sigma^2))^(N/2)*exp(-N/2);
    
    x_removeCaped = x_array;
    %INTENDED TO KEEP ONLY POINTS IN THE TRUST REGION
    x_removeCaped(pointsCaped) = [];
    
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

% CODE RESUMES HERE ONCE EITHER THE LINEAR OR KRIGING MODEL HAS BEEN
% CONSTRUCTED

% THE HYBRID FUNCTION IS PLOTTED HERE ALONG WITH TEST POINTS
plot(testpoints,fhat_h(testpoints),'k','LineWidth',2)
axis([0,inf,-10,25]);
scatter(x_array, yFh_array,200,'o')
scatter(x_array, yFl_array,200,'o')

% DEFINES TRUST REGION AND PLOTS IT
LB = [center-TRS/2]; 
if LB < xMin
    LB = xMin;
end
UB = [center+TRS/2]; 
if UB > xMax
    UB = xMax;
end
plot([LB,UB],[-10,-10],'g','linewidth',2)

% OPTIMIZATION STAGE, FINDS THE MIN LOCAL VALUE OF THE HYBRID FUNCTION
A = []; B = []; Aeq = []; Beq = [];

options = optimset('Algorithm','interior-point','GradObj','off',...
    'GradConstr','off','TolX',1e-6,'Display','iter', 'PlotFcns',...
    @optimplotfval );

% NEW CENTER VALUE OBTAINED IN THE OPTIMIZATION STEP
[center,objectiveValue] = fmincon(@(x)fhat_h(x),center,A,B,Aeq,Beq,LB,...
UB,[],options);

% VECTOR THAT TRACK VALUES OF THE OBJECTIVE HYBRID FUNCTION 
obj_hist = [obj_hist; objectiveValue];

x_array = [x_array; center]; % New evaluation site

% NEW VALUES CALCULATED FOR THE CHEAP AND EXPENSIVE FUNCTIONS AND THEIR
% DERIVATIVES BASED ON THE NEWLY CALCULATED CETNER VALUE.
Fh = fe(center); % new expensive model evaluation

Fl = fc(center); % new cheap model evaluation, should be able to extract this from the optz

dFh_dx = dfe_dx(center); % grad values
dFl_dx = dfc_dx(center); % grad values

yFh_array = [yFh_array; Fh]; yFl_array = [yFl_array; Fl]; 
dFh_array = [dFh_array; dFh_dx]; dFl_array = [dFl_array; dFl_dx];

% NEXT CALCULATE THE RHO VALUE THAT IS USED TO DETERMINE TRUST REGION SIZE
% FOR NEXT STEP
% maybe a differnt way to do this calculation...
rho = ( yFh_array(ii) - yFh_array(ii+1) ) / ( yFh_array(ii) - objectiveValue );

if rho > 0.8
    TRS = 1.05*TRS; % was 2*TRS This creates errors  .7
elseif rho > 0.5
    TRS = TRS;
else
    TRS = 0.5*TRS;
end

% APPEND NEW VALUES FOR TRS AND RHO
TRS_array = [TRS_array;TRS];
rho_array = [rho_array; rho];

% CHECK FOR CONVERGENCE
ConvergenceCheck = (abs(yFh_array(ii) - yFh_array(ii+1))) / abs(yFh_array(ii));
ii = ii + 1;

end
