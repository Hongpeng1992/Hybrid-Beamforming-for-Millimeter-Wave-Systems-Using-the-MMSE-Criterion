function [ W_D,W_RF,MSE ] = MO_AltMinW( H, NRF ,Vn, W_RF )

[Nt, Ns] = size(H);
%W_RF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
[W_RF,MSE] = sig_manif2(H,Vn,W_RF);
M = W_RF'*(H)*H'*W_RF+ Vn*(W_RF)'*W_RF;
W_D = inv(M)*W_RF'*H;
