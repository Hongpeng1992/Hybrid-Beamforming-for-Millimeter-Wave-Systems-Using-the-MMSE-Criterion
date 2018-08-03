function [V_opt] = TWBMSEopt( H, Vn ,N_s)
[N_r,N_t,N_k] = size(H);
V_opt = zeros(N_t,N_s,N_k);
W_opt = zeros(N_r,N_s,N_k);
MSE = zeros(1,N_k);
 for i = 1:N_k
     [ V_opt(:,:,i),W_opt(:,:,i)] = TMSEopt( H(:,:,i), Vn ,N_s);
 end
 MSE = sum(MSE)/N_k;