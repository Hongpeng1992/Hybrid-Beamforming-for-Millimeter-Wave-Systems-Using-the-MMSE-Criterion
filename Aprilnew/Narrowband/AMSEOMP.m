function [V_RF,V_D] = AMSEOMP (N_RF,H,AT,Vn1)
V_RF =[];
[Nt,Ns] = size(H);
Eyy = H*H'+ Vn1*eye(Nt);
C = Eyy^0.5;
Wi = Eyy^(-1)*H;
B = C*Wi;
VRES = B;
for i = 1 : N_RF
    vi = AT'*C'*VRES;
    vi = diag(vi*vi');
    [a,k] = max(vi);
    V_RF =[V_RF AT(:,k)];
    AT(:,k) = [];
    V_D = (V_RF'*C'*C*V_RF)^(-1)*V_RF'*C'*B;
    RES = B-C*V_RF*V_D;
    VRES = RES/norm(RES,'fro');
end
 
    
    
