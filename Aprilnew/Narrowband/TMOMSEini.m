function [t] = TMOMSEini(N_s,N_RF,H,Vn,W_opt)
n = 1;
MSE =zeros(1,11);
[N_r,N_t] = size(H);
% W_D = complex(rand(N_RF,N_s),rand(N_RF,N_s))/10;  %initialization
W_opt = exp( 1i*unifrnd(0,2*pi,N_r,N_s) );
W_RF =exp( 1i*unifrnd(0,2*pi,N_r,N_RF) );
V_RF = exp( 1i*unifrnd(0,2*pi,N_t,N_RF) );
while(n<=3)
    if (n==1)
        H_u = H'*W_opt;
        tw = trace(W_opt'*W_opt);
    else
        H_u = H'*W_RF*W_D;
    end
    Vn1 = Vn*tw;
    if n ==1
        t1 = clock;
        [ V_D,V_RF,MSEV ] = MO_AltMinW( H_u, N_RF ,Vn1, V_RF );
        t2 = clock;
        t = etime(t2,t1);
    end
    
    H_d = H*V_RF*V_D;
    tv = trace(V_RF*(V_D)*V_D'*V_RF');
    Vn2 = Vn*tv;
    [ W_D,W_RF,MSEW ] = MO_AltMinW( H_d, N_RF ,Vn2, W_RF );
    He = W_D'*W_RF'*H*V_RF*V_D;
    tw = trace(W_D'*(W_RF)'*W_RF*W_D);
    MSE(n) = trace(He * He'- He- He'+ eye(N_s)) + Vn2*tw;
    n = n + 1;
end
V_D = V_D/sqrt(tv);
W_D = W_D*sqrt(tv);