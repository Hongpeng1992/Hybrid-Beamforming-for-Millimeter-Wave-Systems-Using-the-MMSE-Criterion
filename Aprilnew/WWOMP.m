function [V_D,V_RF] = WWOMP (H,N_RF,AT,Vn)
V_RF =[];
[Nt, Ns,N_k] = size(H);
for i = 1:N_k
Eyy = H(:,:,i)*H(:,:,i)'+Vn*eye(Nt);
Wi = Eyy^(-1)*H(:,:,i);
C(:,:,i) = Eyy^(0.5);
B(:,:,i) = C(:,:,i)*Wi;
end
VRES = B;
for l = 1 : N_RF
    Vi = [];
    vi = [];
    for m = 1:N_k
        vi(:,:,m) = AT'*C(:,:,m)'*VRES(:,:,m);
        Vi(:,:,m) = vi(:,:,m)*vi(:,:,m)';
    end
    Vi = sum(Vi,3);
    Vi = diag(Vi);
    [a,k] = max(Vi);
    V_RF =[V_RF AT(:,k)];
    AT(:,k) = [];
    V_D = [];
    for i = 1:N_k
        V_D(:,:,i) = (V_RF'*C(:,:,i)'*C(:,:,i)*V_RF)^(-1)*V_RF'*C(:,:,i)'*B(:,:,i);
        RES = B(:,:,i)-C(:,:,i)*V_RF*V_D(:,:,i);
        VRES(:,:,i) = RES/norm(RES,'fro');
    end
end



