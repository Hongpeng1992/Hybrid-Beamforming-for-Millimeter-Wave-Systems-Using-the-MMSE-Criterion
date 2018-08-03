function [gradient] = newGD1(H,Vn,N_RF,W)
[N_r,N_s,N_k] = size(H);
g= zeros(N_r*N_RF,N_k);
for i = 1:N_k
    g(:,i)= nyearGD(H(:,:,i),Vn(i),N_RF,W);
end
gradient = sum(g,2);