function [gradient] = newGD(H,Vn,N_RF,W)
[N_r,N_s] = size(H);
W = reshape(W,N_r,N_RF);
WW = (W'*W)^(-1);
gradient = zeros(N_r,N_RF);
L = (1/Vn*H'*W*WW*W'*H+eye(N_s))^(-1);
A = L*1/Vn*H'*W*WW;
B = H*L;
C = W*WW*W'* B;
for i = 1:N_r
    for j = 1:N_RF
        gradient(i,j)=trace(A(:,j)*B(i,:)+A(:,j)*C(i,:));
    end
end
gradient = -gradient(:);



        