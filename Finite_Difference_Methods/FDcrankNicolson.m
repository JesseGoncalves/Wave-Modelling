function M = FDcrankNicolson(n,r,q)
% outputs matrix for approximating second derivative using finite
% difference methods
% n = number of grid points
% r = time step / (2*grid spacing squared)
% q = designate which matrix (0 -> lhs, =/= 0 -> rhs)
if q == 0
    v1 = (1 + 2*r)*ones(n-2,1);
    v2 = -r*ones(n-3,1);
    M = diag(v1) + diag(v2,1) + diag(v2,-1);
else
    v1 = (1 - 2*r)*ones(n-2,1);
    v2 = r*ones(n-3,1);
    M = diag(v1) + diag(v2,1) + diag(v2,-1);
end
end
