cd('C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\codes')
addpath("C:\compeconfolder\CEtools", "C:\compeconfolder\CEdemos", "C:\dynare\5.3\matlab", "C:\Program Files\Artelys\Knitro 14.0.0\knitromatlab")

clear
dim           = [10,10];
deg                 = 1;
n_shocks      =10;
dynare dz;
[grid, fspace, shocks, prob,N, bounds] = initialize_dz(dim, deg, n_shocks, M_, oo_);

dyn = guess_dynare_dz(M_, oo_, grid);

policy       = dyn.l;
policy(N+1:2*N)   = riskfree_dz(policy, grid, fspace, M_, shocks, prob);

tol         = 1e-8;
maxiter     = 1e5;
obj         = @(pol_next)EE_dz(pol_next,policy, grid, fspace, M_, shocks, prob);
policy_next = knitro_nlneqs(obj, policy(1:N));
policy_next(N+1:2*N) = riskfree_dz(policy_next, grid, fspace, M_, shocks, prob);
%riskfree(policy_next(1:N), grid, fspace, M_, shocks, prob);
%policy_next(2*N+1:3*N)  = coeff_sdf(policy_next,grid,fspace, M_, shocks, prob);
res             = log(policy_next./policy);
resdiff        = policy_next - policy;
lambda =1;
policy      = lambda*policy_next+(1-lambda)*policy;
ii = 1;
norms = [norm(res, inf)];
while norm(res, inf) > tol && ii < maxiter %&& min(policy(N+1:2*N))>0.75 && min(policy(1:N))>0.1
    %if norms(end)<0.01
    %    lambda = 0.01;
    %else
    %    lambda =0.1;
    %end
    obj         = @(pol_next)EE_dz(pol_next,policy, grid, fspace, M_, shocks, prob);
    policy_next = knitro_nlneqs(obj, policy(1:N));
    policy_next(N+1:2*N) = riskfree_dz(policy_next, grid, fspace, M_, shocks, prob);
%riskfree(policy_next(1:N), grid, fspace, M_, shocks, prob);
    %policy_next(2*N+1:3*N)  = coeff_sdf(policy_next,grid,fspace, M_, shocks, prob);
     res             = log(policy_next./policy);
    policy      = lambda*policy_next + (1-lambda)*policy;
    ii              = ii + 1;
    norms(ii+1,:)   = [norm(res,inf)];
    [norm(res,inf)]
end

coeffs_dz.coeffs_unconstr = [funfitxy(fspace,grid, policy(1:N)), funfitxy(fspace,grid, policy(N+1:2*N))];%, ...
    %funfitxy(fspace,grid, policy(2*N+1:3*N))];
save('coeffs_dz.mat', 'coeffs_dz')

load('coeffs_forward.mat')
coeffs = coeffs_forward.coeffs_02;
theta_mean = log(0.2);
M_.params(end) = theta_mean;
G =20;
xx = [linspace(bounds(1,1),bounds(1,2),G)', linspace(bounds(2,1),bounds(2,2),G)', linspace(bounds(3,1),bounds(3,2),G)'];
grid_new = gridmake(xx(:,1), xx(:,2), xx(:,3));
N = size(grid_new,1);
policy       = funeval(coeffs(:,1), fspace, grid_new);
policy(N+1:2*N)   = funeval(coeffs(:,2), fspace, grid_new);
policy(2*N+1:3*N)   = funeval(coeffs(:,3), fspace, grid_new);
err = EulerErr(policy, grid_new, fspace, M_, shocks, prob);