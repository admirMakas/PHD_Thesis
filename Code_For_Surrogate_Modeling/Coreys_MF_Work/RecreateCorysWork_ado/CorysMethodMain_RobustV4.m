%clear all
close all
clc
warning off

%% inputs 
TYPE = 'DAN'; % 'COR'
%Bounds = [0;1.1];
Bounds = [-2, -2; 2, 2];
%Center = [0.55];
Center = [-2, -2];


gridSize = 2000; % 'Display','iter', interior-point
options = optimset('Algorithm','sqp','GradObj','on',...
    'GradConstr','off','TolX',1e-12, 'PlotFcns',...
    @optimplotfval );

[ Hist ] = nDimensionalMF( TYPE, Center, Bounds, options, 'YES',gridSize);

% [center,objectiveValue,~,OUTPUTCON] = fmincon(...
%             @(x)HighFidelity( x ),Center,[],[],[], [],Bounds(1,:), Bounds(2,:),[],options);
%         

% This is the observed solution
% load '1D_solution.mat'       