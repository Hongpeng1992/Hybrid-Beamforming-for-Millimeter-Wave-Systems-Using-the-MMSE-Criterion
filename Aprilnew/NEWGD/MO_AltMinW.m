function [ W_D,W_RF,MSE ] = MO_AltMinW( H, N_RF ,Vn, W_RF )

[Nt, Ns,N_k] = size(H);
%W_RF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
[W_RF,MSE] = sig_manif2(H,Vn,W_RF);
W_D = zeros(N_RF,Ns,N_k);
for i = 1:N_k
M = W_RF'*(H(:,:,i))*H(:,:,i)'*W_RF+ Vn(i)*(W_RF)'*W_RF;
W_D(:,:,i) = inv(M)*W_RF'*H(:,:,i);
end
MSE = MSE/N_k;
