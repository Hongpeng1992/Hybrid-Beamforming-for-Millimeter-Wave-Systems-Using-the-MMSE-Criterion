function [gradient] = newGD(H,Vn,N_RF,W)
[N_r,N_s] = size(H);
W = reshape(W,N_r,N_RF);
WW = (W'*W)^(-1);
gradient = zeros(N_r,N_RF);
A = (1/Vn*H'*W*WW*W'*H+eye(N_s))^(-2);
B = H'*W*WW;
C = B';
M = 1/Vn*W*C*A*B;
N = 1/Vn*H*A*B;
gradient = M-N;
gradient = gradient(:);
        