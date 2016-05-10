%This code will run example 2.6 from book. Effective pages are
%47-48

%This example pertains to 

%clear workspace
clear all
clc

% Make ModelInfo visible to all functions
global ModelInfo

% Sampling plan
ModelInfo.X=bestlh(10,2,50,25);

% Compute objective function values – in this case using
% the dome.m test function
for i=1:size(ModelInfo.X,1)
ModelInfo.y(i)=dome(ModelInfo.X(i,:));
end
ModelInfo.y=ModelInfo.y';

% Basis function type:
ModelInfo.Code=3;

% Estimate model parameters
rbf

% Plot the predictor
x=(0:0.025:1);
for i=1:length(x)
for j=1:length(x)
    M(j,i)=predrbf([x(i) x(j)]);
end
end
contour(x,x,M)

hold on
for i=1:length(x)
for j=1:length(x)
    N(j,i)=dome([x(i) x(j)]);
end
end
contour(x,x,N)

figure
contour3(x,x,M,30,'-k', 'LineWidth', 2)
hold on
contour3(x,x,N,30,'--r', 'LineWidth', 2)
