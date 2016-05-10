clear all
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
Xe = [0, .3,0.4, 0.6, 1]'; Re = fe(Xe); % These points most coincide with Xc
Xc = Xe; Rc = fc(Xc); dRe = dfe_dx(Xe);

%% My code
[ model ] = gekFit( Xe,Xe, Re, dRe);
[yhat_mine, MSE] = gekPred( model,testpoints );

%% Other code
options.debug = 2;
[yhat_not, varxout, gradout, vargradout, report] = ...
    gek(Xe,Re,zeros(size(Re)),dRe,zeros(size(dRe)),testpoints',options);



    %% Plot for each run
figure, hold on
scatter(Xe, Re)
plot(testpoints,fe(testpoints),'k','LineWidth',2)
plot(testpoints,yhat_mine,'b','LineWidth',2)
plot(testpoints,yhat_not,'r','LineWidth',2)
legend('Data Points','Y_{Expensive}','My Gek','Other Gek') 

