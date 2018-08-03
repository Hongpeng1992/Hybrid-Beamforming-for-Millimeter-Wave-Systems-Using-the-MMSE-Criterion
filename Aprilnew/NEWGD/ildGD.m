function [gradient] = newGD(H,Vn,N_RF,W)
t1 = cputime;
for m = 1:10000
[N_r,N_s] = size(H);
W = reshape(W,N_r,N_RF);
WW = (W'*W)^(-1);
gradient = zeros(N_r,N_RF);
A = (1/Vn*H'*W*WW*W'*H+eye(N_s))^(-1);
for i = 1:N_r
    for j = 1:N_RF
        M = zeros(N_RF,N_r);
        M(j,i) = 1;
        B = 1/Vn*H'*W*(WW*M*W*WW*W'+WW*M)*H;
        gradient(i,j)=trace(A*B*A);
    end
end
gradient = -gradient(:);
end
t2 = cputime;
t = t2-t1;
        