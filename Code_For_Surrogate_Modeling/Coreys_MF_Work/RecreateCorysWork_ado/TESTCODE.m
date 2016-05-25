clear all
close all
clc

thetas = [.5,.5,.5]
x = [0, 0, 0; 
    2 , 2, 2]
S = [1, 1, 1]



xC = repmat(x,length(thetas)+1,1)
sC = repmat(S,length(xC),1)
thetasC = repmat(thetas,length(xC),1)
r = exp(-1*sum(thetasC.*abs(xC-sC).^2,2))
r = kron(r,ones(1,3) )



    

