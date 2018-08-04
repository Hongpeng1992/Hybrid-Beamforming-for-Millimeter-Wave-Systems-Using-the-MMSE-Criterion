function [V_D,V_RF,W_D,W_RF] = JZPE(N_RF,V_opt,W_opt)
  %% transmit 
     a = norm(W_opt,'fro');
     [ V_RF,V_D ] = PE_AltMin( V_opt, N_RF );
     [ W_RF,W_D ] = PE_AltMin( W_opt, N_RF );
     V_D = V_D / norm(V_RF*V_D, 'fro') ; 
     W_D = W_D / norm(W_RF*W_D, 'fro') * a; 