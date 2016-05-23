clear all
close all
clc
warning off

%% inputs 
TYPE = 'DAN'; % 'COR'
Bounds = [0;1.1];
gridSize = 1000; % 'Display','iter',
options = optimset('Algorithm','interior-point','GradObj','on',...
    'GradConstr','off','TolX',1e-6, 'PlotFcns',...
    @optimplotfval );

[ Hist ] = nDimensionalMF( TYPE, [0.55], Bounds, options, 'YES',gridSize);

% [center,objectiveValue,~,OUTPUTCON] = fmincon(...
%             @(x)HighFidelity( x ),.3,[],[],[], [],0, 1.1,[],options);
%         

% This is the observed solution
% load '1D_solution.mat'       