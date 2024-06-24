cd('C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\codes')
addpath("C:\compeconfolder\CEtools", "C:\compeconfolder\CEdemos", ...
    "C:\dynare\5.3\matlab", "C:\Program Files\Artelys\Knitro 14.0.0\knitromatlab",...
    "C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\coeffs",...
    "C:\Users\giaco\Desktop\phd\myfunc")

clear
dim           = [10,10];
deg                 = 1;
n_shocks      =5;
dynare dz;
[grid, fspace, shocks, prob,N, bounds] = initialize_dz(dim, deg, n_shocks, M_, oo_);
inve =0 ;
dyn = guess_dynare_dz(M_, oo_, grid);

policy       = log(dyn.l);
policy(N+1:2*N)   = riskfree_dz(policy, grid, fspace, M_, shocks, prob, inve);
%M_.params(2) = 4;
tol         = 1e-5;
maxiter     = 1e5;
obj         = @(pol_next)EE_dz_unconstr(pol_next,policy, grid, fspace, M_, shocks, prob, inve);
policy_next = knitro_nlneqs(obj, policy(1:N));
policy_next(N+1:2*N) = riskfree_dz(policy_next, grid, fspace, M_, shocks, prob, inve);
res             = log(policy_next./policy);
resdiff        = policy_next - policy;
lambda = 0.1;
policy      = lambda*policy_next+(1-lambda)*policy;
ii = 1;
norms = [norm(res, inf)];
while norm(res, inf) > tol && ii < maxiter
    obj         = @(pol_next)EE_dz_unconstr(pol_next,policy, grid, fspace, M_, shocks, prob, inve);
    policy_next = knitro_nlneqs(obj, policy(1:N));
    policy_next(N+1:2*N) = riskfree_dz(policy_next, grid, fspace, M_, shocks, prob, inve);
    res             = log(policy_next./policy);
    policy      = lambda*policy_next + (1-lambda)*policy;
    ii              = ii + 1;
    norms(ii+1,:)   = [norm(res,inf)];
    norm(res,inf)
end

coeffs_dz_.coeffs_unconstr = [funfitxy(fspace,grid, policy(1:N)), funfitxy(fspace,grid, policy(N+1:2*N))];
save('coeffs_dz_cal2.mat', 'coeffs_dz_cal2')