function [W_RF,cost] = sig_manif2(H,Vn,W_RF)
[Nt, NRF] = size(W_RF);
manifold = complexcirclefactory(Nt*NRF);
problem.M = manifold;
[x,cost,info,reason,options] = conjugategradient1(problem,W_RF(:),H,Vn,NRF);
W_RF = reshape(x,Nt,NRF);
end