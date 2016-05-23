function [ Hist ] = EvaluateFidelities( center, Hist )


[ Fh, dFh_dx ] = HighFidelity( center );
[ Fl, dFl_dx ] = LowFidelity( center );

Hist.yFh = [Hist.yFh; Fh]; Hist.yFl = [Hist.yFl; Fl]; 
Hist.dFh_dx = [Hist.dFh_dx; dFh_dx]; Hist.dFl_dx = [Hist.dFl_dx; dFl_dx];


end

