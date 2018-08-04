function [V_D,V_RF,W_D,W_RF,MSE,n] = AOMPini(Ns,N_RF,H,Vn,W_opt,AT,AR)
 n = 1;   %number of iteration
 MSE =zeros(1,11);
 W_D = W_opt;  %initialization
 W_RF = 1;
 tw = trace(W_D'*(W_RF)'*W_RF*W_D);
 [~,N_t] = size(H);
 %V_RF = exp( 1i*unifrnd(0,2*pi,N_t,N_RF) );
 tv=1;
 while(n<3 || (MSE(n-2)-MSE(n-1))>1e-4 &&n<=10)
     H_u = H'*W_RF*W_D;          %effective downlink channel   
     Vn1 = tw*Vn;
     [V_RF,V_D] = MSEOMP (N_RF,H_u,AT,Vn1);
     tv = trace(V_RF*(V_D)*V_D'*V_RF');     
     %%UEside
     V_D = V_D/sqrt(tv);
     H_d = H*V_RF*V_D;
     [W_RF,W_D] = AMSEOMP (N_RF,H_d,AR,Vn); %the same formulation
     beta = tv^-0.5;
     tw = trace(W_D'*(W_RF)'*W_RF*W_D);     
     Hequal= W_D'*W_RF'*H*V_RF*V_D;
     MSE(n)=  trace(beta^-2*Hequal*Hequal'-beta^-1*Hequal-beta^-1*Hequal'+Vn*beta^-2*W_D'*W_RF'*W_RF*W_D+eye(Ns));
     n = n + 1;     
 end
       
       
     %  W_D = W_D*sqrt(tv);