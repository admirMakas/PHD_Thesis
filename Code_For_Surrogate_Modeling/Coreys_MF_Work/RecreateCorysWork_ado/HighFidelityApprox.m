function [ fhat_h, dfhat_h_dx ] = HighFidelityApprox( x, modelAdd, modelMult, w1, w2 )
    [ADD, ~, dADD_dx] = gekPred( modelAdd,x );
    [MULT, ~, dMULT_dx] = gekPred( modelMult,x );
    [ Fl, dFl_dx ] = LowFidelity( x );

fhat_h =  w1.*(Fl + ADD ) + w2.*( Fl.*MULT );

if nargout == 1
    dfhat_h_dx = [];
    
else
    
dfhat_h_dx = w1.*dFl_dx + w1.*dADD_dx + w2.*(Fl.*dMULT_dx + dFl_dx.*MULT);

end


end

