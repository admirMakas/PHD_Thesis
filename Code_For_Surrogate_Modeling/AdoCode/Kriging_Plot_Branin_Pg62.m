Xplot=0:1/20:1;

for i=1:21
	for j=1:21
        BraninPred(j,i)=pred([Xplot(i) Xplot(j)]);
        BraninTrue(j,i)=branin([Xplot(i) Xplot(j)]);
    end
end

%Contour plot
figure
contour(Xplot, Xplot, BraninPred, 35)
hold on
contour(Xplot, Xplot, BraninTrue, 35, '-r')
hold on

scatter(ModelInfo.X(:,1), ModelInfo.X(:,2))