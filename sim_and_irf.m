cd('C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\codes')
addpath("C:\compeconfolder\CEtools", "C:\compeconfolder\CEdemos", ...
    "C:\dynare\5.3\matlab", "C:\Program Files\Artelys\Knitro 14.0.0\knitromatlab",...
    "C:\Users\giaco\Desktop\phd\Current Projects\Uncertaieps_thetanty\coeffs01",...
     "C:\Users\giaco\Desktop\phd\myfunc")

clear
dim           = [12,10];
deg           = 1;
n_shocks      = 5;
dynare dz;
[grid, fspace, shocks, prob,N, bounds, pars] = initialize_dz(dim, deg, n_shocks, M_, oo_);
inve = 0; %0 = convex cost , 1 = concave prod func

load('coeffs_dz_rightreturn.mat')
coeffs01 = coeffs_dz_.coeffs_381_01;
coeffs1 = coeffs_dz_.coeffs_381_1;
X0 = [mean(grid(:,1)), 0];
policy       = funeval(coeffs01(:,1), fspace, X0);
policy(2)   = funeval(coeffs01(:,2), fspace, X0);
policy(3)   = funeval(coeffs01(:,3), fspace, X0);
pars(end-1) = 0.381;
pars(end) = 0.01;
T1 = 100;
T2 = 100;
nn = 500;
simul_shocks_z = randn(T1, nn);
simul_shocks_d = randn(T1, nn);
ergo_ss = zeros(5,nn);
init_state=zeros(nn,2);
for n = 1:nn
    X = [X0(1)+pars(end)*simul_shocks_d(1,n), X0(2)+sd_z*simul_shocks_z(1,n)];
    for t = 1:T1
        pf = policy_functions_for_sim(coeffs01, X, pars, fspace,  inve);
        X = [pf.GG_prime(1)+pars(end)*simul_shocks_d(1,n), pf.GG_prime(2)+sd_z*simul_shocks_z(1,n)];
    end
    ergo_ss(:,n) = [log(pf.l);pf.R;max(pf.psi,1e-6); pf.y; pf.Q];
    init_state(n,:) = pf.GG_prime;
end

unc = 0.01;
pars(end) = unc;
irf0 = zeros(T2,5,nn);
for n = 1:nn
    X = init_state(n,:);
    for t = 1:T2
        pf = policy_functions_for_sim(coeffs01, X, pars, fspace,  inve);
        X = [pf.GG_prime(1), pf.GG_prime(2)];
        irf0(t,:,n) = [pf.l, pf.R, pf.psi ,pf.y, pf.Q];
    end
end

unc = 0.1;
pars(end) = unc;
irf1 = zeros(T2,5,nn);
for n = 1:nn
    X = init_state(n,:);
    for t = 1:T2
        pf = policy_functions_for_sim(coeffs1, X, pars, fspace,  inve);
        X = [pf.GG_prime(1), pf.GG_prime(2)];
        irf1(t,:,n) = [pf.l, pf.R, pf.psi ,pf.y, pf.Q];
        %pars(end) = 0.5*pars(end) + 0.5*0.01;
    end
end

irf0_mean = sum(irf0, 3)/nn;
irf1_mean = sum(irf1, 3)/nn;
irf1_0 = log(irf1_mean) - log(irf0_mean);
irf_lab = real(irf1_0(:,5));
plot(irf_lab)