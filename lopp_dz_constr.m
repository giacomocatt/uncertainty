cd('C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\codes')
addpath("C:\compeconfolder\CEtools", "C:\compeconfolder\CEdemos", ...
    "C:\dynare\5.3\matlab", "C:\Program Files\Artelys\Knitro 14.0.0\knitromatlab",...
    "C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\coeffs",...
     "C:\Users\giaco\Desktop\phd\myfunc")

clear
dim           = [10,10];
deg           = 1;
n_shocks      = 5;
dynare dz_constr;
[grid, fspace, shocks, prob,N, bounds] = initialize_dz(dim, deg, n_shocks, M_, oo_);
inve = 0; %0 = convex cost , 1 = concave prod func

load('coeffs_dz_lin_inv.mat')
coeffs = coeffs_dz_.coeffs_03;
policy       = funeval(coeffs(:,1), fspace, grid);
policy(N+1:2*N)   = funeval(coeffs(:,2), fspace, grid);
policy(2*N+1:3*N)   = funeval(coeffs(:,3), fspace, grid);%beta*(1-sigma)/(1-beta*sigma)*ones(N,1);
unc = 0.01;
M_.params(end) = unc;
%for jj =1:10
theta_mean = 0.0325;
M_.params(end-1) = theta_mean;
%M_.params(2) = 1;
tol         = 1e-5;
maxiter     = 1e5;
obj         = @(pol_next)EE_dz(pol_next,policy, grid, fspace, M_, shocks, prob, inve);
policy_next = knitro_nlneqs(obj, policy(1:N));
policy_next(N+1:3*N) = R_psi_dz(policy, policy_next,grid,fspace, M_, shocks, prob, inve);
res             = log(policy_next./policy);
resdiff        = policy_next - policy;
lambda = 0.1;
policy      = lambda*policy_next+(1-lambda)*policy;
ii = 1;
norms = [norm(res, inf)];
while norm(res, inf) > tol && ii < maxiter %&& min(policy(N+1:2*N))>0.75 && min(policy(1:N))>0.1
    obj         = @(pol_next)EE_dz(pol_next,policy, grid, fspace, M_, shocks, prob, inve);
    policy_next = knitro_nlneqs(obj, policy(1:N));
    policy_next(N+1:3*N) = R_psi_dz(policy, policy_next,grid,fspace, M_, shocks, prob, inve);
    res             = log(policy_next./policy);
    policy      = lambda*policy_next + (1-lambda)*policy;
    ii              = ii + 1;
    norms(ii+1,:)   = [norm(res,inf)];
    norm(res,inf), M_.params(end-1)
end
store = policy;
end
policy = store;
coeffs_dz_.coeffs_03 = [funfitxy(fspace,grid, policy(1:N)), funfitxy(fspace,grid, policy(N+1:2*N)), ...
    funfitxy(fspace,grid, policy(2*N+1:3*N))];
save('coeffs_dz_convex_inv.mat', 'coeffs_dz_')

dim2 = [20,10];
[grid1, fspace1, Phi] = Grid_spli(dim2,bounds,deg);
pol_func = policy_functions_dz(coeffs, grid, M_, fspace,fspace, shocks, prob, inve);