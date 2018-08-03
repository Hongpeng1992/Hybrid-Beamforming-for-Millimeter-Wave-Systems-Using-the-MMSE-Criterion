function [V_D,V_RF,W_D,W_RF,MSE1,MSE] = WBJUNZhang(N_RF,V_opt,W_opt,Vn,H)
[N_r,N_s,N_k] = size(W_opt);
[V_RF,V_D] = MO_AltMinWB( V_opt, N_RF );
W_opt = W_opt*20;  %for better performance
[W_RF,W_D] = MO_AltMinWB( W_opt, N_RF );
W_D = W_D/20;
Rnn = Vn *eye(N_r) ;
Hequal = zeros(N_s,N_s,N_k);
MSE = zeros(1,N_k);
tv = zeros(1,N_k);
for i = 1:N_k
tv(i) = trace(V_RF*(V_D(:,:,i))*V_D(:,:,i)'*V_RF');
V_D(:,:,i) = V_D(:,:,i)/sqrt(tv(i));
Hequal(:,:,i)=W_D(:,:,i)'*W_RF'*H(:,:,i)* V_RF*V_D(:,:,i);
MSE(i) = trace((Hequal(:,:,i) - eye(N_s)) * (Hequal(:,:,i) - eye(N_s))' + W_D(:,:,i)'*W_RF'*Rnn*W_RF*W_D(:,:,i));
end
MSE1 = sum(MSE)/N_k;