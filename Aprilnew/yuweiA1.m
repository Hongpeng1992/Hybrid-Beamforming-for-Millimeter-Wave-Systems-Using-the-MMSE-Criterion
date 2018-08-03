function V_RF  = yuweiA1(Vn,H,NRF,Nt)
V_RF = ones(Nt,NRF);
Nk = size(H,3);
for k = 1:Nk
    F(:,:,k) = H(:,:,k)'*H(:,:,k);
end
F = sum(F,3)/Nk;
g = 1/NRF/Nt;
a = g/Vn;

for Nloop = 1:10
    for j = 1:NRF
        VRF = V_RF;
        VRF(:,j)=[];
        C = eye(NRF-1)+a*VRF'*F*VRF;
        G = a*F-a^2*F*VRF*C^(-1)*VRF'*F;
        for i = 1:Nt
            for l = 1:Nt
                if i~=l
                    x(l)=G(i,l)*V_RF(l,j);
                end
            end
            n = sum(x);
            if n ==0
                V_RF(i,j)=1;
            else
                V_RF(i,j)=n/abs(n);
            end
        end
    end
end