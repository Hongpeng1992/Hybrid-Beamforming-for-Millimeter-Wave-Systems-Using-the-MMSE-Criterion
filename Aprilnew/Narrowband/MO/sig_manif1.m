function [W_RF, cost] = sig_manif1(H, W_D, NRF , tv,Vn,W_RF)
[Nt, NRF] = size(W_RF);

manifold = complexcirclefactory(Nt*NRF);
problem.M = manifold;

% problem.cost  = @(x) norm( Fopt - reshape(x,Nt,NRF) * FBB,'fro')^2;
% problem.egrad = @(x) -2 * kron(conj(FBB), eye(Nt)) * (Fopt(:) - kron(FBB.', eye(Nt)) * x);
h = H(:);
A = kron(W_D.', eye(Nt));
B = kron(W_D.',H');

problem.cost  = @(x) x'*B'*B*x-h'*A*x-x'*A'*h+Vn*tv*x'*A'*A*x+NRF;
problem.egrad = @(x) B'*B*x-A'*h+Vn*tv*A'*A*x;  

% checkgradient(problem);
warning('off', 'manopt:getHessian:approx');

[x,cost,info,options] = conjugategradient(problem,W_RF(:));
% [x,cost,info,options] = trustregions(problem, FRF(:));
% info.iter
W_RF = reshape(x,Nt,NRF);

end