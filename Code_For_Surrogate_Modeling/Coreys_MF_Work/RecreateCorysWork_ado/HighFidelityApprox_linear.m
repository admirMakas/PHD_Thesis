function [ fhat_h, dfhat_h_dx ] = HighFidelityApprox_linear( x, Hist, intercepts, pointsCaped, w1, w2 )

interceptAdd = intercepts(1,:);
interceptMult = intercepts(2,:);

[ADD,dADD_dx] = LinearModel_dan( x, Hist.centers(end,:),  Hist.yAdd(pointsCaped), ...
    Hist.dAdd_dx(pointsCaped), interceptAdd);

[MULT,dMULT_dx] = LinearModel_dan( x, Hist.centers(end,:), Hist.yMult(pointsCaped), ...
    Hist.dMult_dx(pointsCaped),interceptMult);

[ Fl, dFl_dx ] = LowFidelity( x );
fhat_h =  w1.*(Fl + ADD ) + w2.*( Fl.*MULT );
dfhat_h_dx = w1.*dFl_dx + w1.*dADD_dx + w2*(Fl.*dMULT_dx + dFl_dx.*MULT);




end

