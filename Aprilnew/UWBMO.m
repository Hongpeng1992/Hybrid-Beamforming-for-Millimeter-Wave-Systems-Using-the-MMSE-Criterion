function [V_D,V_RF,W_D,W_RF,MSE1] = UWBMO(N_s,N_RF,H,Vn,W_opt)
[N_r,N_t,N_k] = size(H);
%Ëæ»ú³õÊ¼»¯
V_RF = exp( 1i*unifrnd(0,2*pi,N_t,N_RF) );
W_RF = exp( 1i*unifrnd(0,2*pi,N_r,N_RF) );
W_opt = exp( 1i*unifrnd(0,2*pi,N_r,N_s,N_k) ); 
H_u = zeros(N_t,N_s,N_k);
Vn1 = zeros(1,N_k);
H_d = zeros(N_r,N_s,N_k);
Vn2 = zeros(1,N_k);
tv = zeros(1,N_k);
Rnn = Vn *eye(N_r) ;
n = 1;
Hequal = zeros(N_s,N_s,N_k);
MSE = zeros(1,N_k);
MSE1 = zeros(1,10);
while(n<3 || (MSE1(n-2)-MSE1(n-1))>1e-4 &&n<=10 )
    if n ==1
        for i = 1:N_k
            H_u(:,:,i) = H(:,:,i)'*W_opt(:,:,i);
            Vn1(i) = Vn*trace(W_opt(:,:,i)'*W_opt(:,:,i));
        end
    else
        for i = 1:N_k
            H_u(:,:,i) = H(:,:,i)'*W_RF*W_D(:,:,i);
            Vn1(i) = Vn*trace(W_D(:,:,i)'*W_RF'*W_RF*W_D(:,:,i));
        end
    end
    [ V_D,V_RF ] = MO_AltMinW( H_u, N_RF ,Vn1, V_RF );
    for i = 1:N_k
        H_d(:,:,i) = H(:,:,i)*V_RF*V_D(:,:,i);
        tv(i) = trace(V_RF*(V_D(:,:,i))*V_D(:,:,i)'*V_RF');
        Vn2(i) = Vn * tv(i);
    end
    [ W_D,W_RF] = MO_AltMinW( H_d, N_RF ,Vn2, W_RF );
    for i = 1:N_k
        V_D(:,:,i) = V_D(:,:,i)/sqrt(tv(i));
        W_D(:,:,i) = W_D(:,:,i)*sqrt(tv(i));
        Hequal(:,:,i)=W_D(:,:,i)'*W_RF'*H(:,:,i)* V_RF*V_D(:,:,i);
        MSE(i) = trace((Hequal(:,:,i) - eye(N_s)) * (Hequal(:,:,i) - eye(N_s))' + W_D(:,:,i)'*W_RF'*Rnn*W_RF*W_D(:,:,i));
    end
    MSE1(n) = sum(MSE)/N_k;
    n = n + 1;
end