function [ Hist ] = AdditiveAndMultiplicative...
    ( Itr, Hist )

Add = Hist.yFh(Itr) - Hist.yFl(Itr);
dAdd_dx = Hist.dFh_dx(Itr) - Hist.dFl_dx(Itr);
Mult = Hist.yFh(Itr)/Hist.yFl(Itr);
dMult_dx = (Hist.yFl(Itr).*Hist.dFh_dx(Itr) - ...
    Hist.yFh(Itr).*Hist.dFl_dx(Itr))./Hist.yFl(Itr).^2;


Hist.yAdd = [Hist.yAdd; Add]; Hist.yMult = [Hist.yMult; Mult]; 
Hist.dAdd_dx = [Hist.dAdd_dx; dAdd_dx]; Hist.dMult_dx = [Hist.dMult_dx;...
    dMult_dx];


end


%% old code
% function [ Add, Mult, dAdd_dx, dMult_dx ] = AdditiveAndMultiplicative...
%     ( Fh, Fl, dFh_dx, dFl_dx )
% 
% Add = Fh - Fl;
% dAdd_dx = dFh_dx - dFl_dx;
% Mult = Fh./Fl;
% dMult_dx = (Fl.*dFl_dx - Fh.*dFh_dx)./Fl.^2;
% 
% end
