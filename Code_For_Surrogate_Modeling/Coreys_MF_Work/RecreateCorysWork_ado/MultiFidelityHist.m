function [ Hist ] = MultiFidelityHist( )

Hist = struct('centers',[], 'yAdd',[], 'yMult',[], ...
    'dAdd_dx',[], 'dMult_dx',[], 'yFh',[], 'yFl', [], 'dFh_dx',[],...
    'dFl_dx', [], 'TRS', [], 'rho', [] ,'obj' , [], ...
    'w1', [], 'w2', [],'optzLcount',0,'AddModels',[],'MultModels',[],...
    'HighModels',[]);

end

