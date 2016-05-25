close all
clear all
clc

ROS = @(x) 100*(x(:,2)-x(:,1).^2).^2+(1-x(:,1)).^2;

dROS_dx = @(x) [-400*(x(:,2)-x(:,1).^2).*x(:,1)-2.*(1-x(:,1)) , 200.*(x(:,2)-x(:,1).^2)];



Samples = lhsdesign(5,2);
a= -2; b = 2;

Samples(:,1) = a + (b-a).*Samples(:,1);
Samples(:,2) = a + (b-a).*Samples(:,2);

Y = ROS(Samples);

DY = dROS_dx(Samples);
[ model ] = gekFit( Samples, Samples, Y, DY);

gridSize = 30;
testpoints = gridsamp([-2,-2;2,2], gridSize);
Ytrue = ROS(testpoints);

X1 = reshape(testpoints(:,1),gridSize,gridSize);
X2 = reshape(testpoints(:,2),gridSize,gridSize);

ResponseTrue = reshape(Ytrue, size(X1));
[yhat,MSE,dyhat_dx] = gekPred( model,testpoints );
YYY = reshape(yhat,gridSize,gridSize);

figure, hold on
contour(X1,X2,YYY,50)
scatter(Samples(:,1),Samples(:,2),50,'Fill','k')

figure, hold on
contour(X1,X2,ResponseTrue,50)
scatter(Samples(:,1),Samples(:,2),50,'Fill','k')

[a,b,c] = gekPred( model,Samples )
