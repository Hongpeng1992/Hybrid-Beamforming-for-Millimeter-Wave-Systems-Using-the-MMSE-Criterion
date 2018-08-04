function FL = LASSO_soft( FBD,FBB )
[Nt,Ns,K,F] = size(FBD);
Fopt = reshape(FBD,[Nt,Ns*K*F]);

temp = Fopt*FBB';
FRF = temp - exp(1i*angle(temp)).*max(0,abs(temp)-2);

FL = FRF*FBB;
if(norm(FL,'fro')^2>K*Ns*F)
    FL = sqrt(K*F*Ns)*FL/norm(FL,'fro');
end
FL = reshape(FL,[Nt,Ns,K,F]);
end