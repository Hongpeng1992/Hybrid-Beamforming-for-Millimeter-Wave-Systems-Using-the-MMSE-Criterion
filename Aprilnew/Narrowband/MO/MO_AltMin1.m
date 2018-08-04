function [ W_D,W_RF,MSE ] = MO_AltMin1( H, NRF ,tv,Vn,W_RF )

[Nt, Ns] = size(H);
y = [];
W_RF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
while(isempty(y) || (y(1)-y(2))>1e-5)
    M = W_RF'*(H)*H'*W_RF+ Vn*tv*(W_RF)'*W_RF;
    W_D = inv(M)*W_RF'*H;
    Hequal = W_D'*W_RF'*H;
    y(1) = trace((Hequal - eye(NRF)) * (Hequal - eye(NRF))' + W_D'*W_RF'* Vn *eye(Nt)*W_RF*W_D);
    %y(1) = norm(Fopt - FRF * FBB,'fro')^2;
   
    [W_RF , y(2)] = sig_manif1(H, W_D, NRF , tv,Vn,W_RF);
end
    MSE = y(2);
end