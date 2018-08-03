function [V_D,V_RF] = WOMP (H,N_RF,AT,Vn)
V_RF =[];
[Nt, Ns,N_k] = size(H);
for i = 1:N_k
    VRES(:,:,i)  = eye(Ns);
end
for l = 1 : N_RF
    Vi = [];
    vi = [];
    for m = 1:N_k
        vi(:,:,m) = AT'*H(:,:,m)*VRES(:,:,m);
        Vi(:,:,m) = vi(:,:,m)*vi(:,:,m)';
    end
    Vi = sum(Vi,3);
    Vi = diag(Vi);
    [a,k] = max(Vi);
    V_RF =[V_RF AT(:,k)];
    AT(:,k) = [];
    V_D = [];
    for i = 1:N_k
        M = V_RF'*(H(:,:,i))*H(:,:,i)'*V_RF+ Vn(i)*(V_RF)'*V_RF;
        V_D(:,:,i) = inv(M)*V_RF'*H(:,:,i);
        RES = eye(Ns)-H(:,:,i)'*V_RF*V_D(:,:,i);
        VRES(:,:,i) = RES/norm(RES,'fro');
    end
end



