global ModelInfo

% Number of variables
k=2;
% Number of sample points
n=20;

% Create sampling plan
ModelInfo.X=bestlh(n,k,50,20);

% Calculate observed data
for i=1:n
    ModelInfo.y(i,1)=branin(ModelInfo.X(i,:));
end

% Set upper and lower bounds for search of log theta
UpperTheta=ones(1,k).*2;
LowerTheta=ones(1,k).*-3;

% Run GA search of likelihood
[ModelInfo.Theta,MinNegLnLikelihood]=...
ga(@likelihood,k,[],[],[],[], LowerTheta,UpperTheta);

% Put Cholesky factorization of Psi, into ModelInfo
[NegLnLike,ModelInfo.Psi,ModelInfo.U]=likelihood(ModelInfo.Theta);