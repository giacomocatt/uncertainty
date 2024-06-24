function [pol_func] = policy_functions_dz(coeffs, grid, M_, fspace, fspace1, shocks, prob, inve)
beta        = M_.params(1);       
gamma       = M_.params(2);  
chi         = M_.params(3);
delta       = M_.params(4);  
alpha       = M_.params(5);       
a1          = M_.params(6);                  
csi         = M_.params(7);   
zeta        = M_.params(8); 
f           = M_.params(9);
omega       = M_.params(10);
sigma       = M_.params(11);
z_mean      = M_.params(12);
sd_z        = M_.params(13); 
rho_z       = M_.params(16);

theta_mean  = M_.params(end-1);
unc         = M_.params(end);
N          = size(grid,1);
policy(1:N) = funeval(coeffs(:,1), fspace, grid);
policy(N+1:2*N)   = funeval(coeffs(:,2), fspace, grid);
policy(2*N+1:3*N)   = funeval(coeffs(:,3), fspace, grid);

l        =  exp(policy(1:N))';

R           =  policy(N+1:2*N)';
psi     = policy(2*N+1:3*N)';

Q = (chi*l - (1-f)*(1-alpha)*exp(grid(:,2)).*l.^(1-alpha))/(zeta*f);
y       = exp(grid(:,2)).*(l).^(1-alpha);
       
i       = (Q/a1).^(1/csi);

if inve == 0
    gk = 1-delta + i;
    cost =  a1*i.^(1+csi)/(1+csi)-csi*delta/(1+csi);
else
    gk = 1-delta + a1^(-1)*i.^(1-csi)/(1-csi)-csi*delta/(1-csi);
    cost = i;
end

c       =  max(y-cost,0.01);

Rk      =  gk.*Q+alpha*y -cost ;

n_net  = max(Rk - exp(grid(:,1)),0.1);

n_prime = (sigma*(n_net) + omega*Q)./(gk);

d_prime = R.*(Q - n_prime);

u       = c - chi*l;

constr = theta_mean*Q./(n_prime);

pol_func.l = l;
pol_func.R = R;
pol_func.psi = psi;
pol_func.Q = Q;
pol_func.y = y;
pol_func.c = c;
pol_func.i = i;
pol_func.Rk = Rk./Q;
pol_func.n_net = n_net;
pol_func.n_prime = n_prime;
pol_func.d_prime = d_prime;
pol_func.u = u;

muc     = u.^(-gamma);
sdf     = psi.*muc;
E_Rk    = Rk;

grid_prime   = [log(d_prime), rho_z*grid(:,2)+(1-rho_z)*z_mean];
E = expect_coeff_dz(grid_prime,unc, sd_z, shocks, prob, fspace1);
pol_func.E = E;

gamma_lam   = funfitxy(fspace1, grid, muc);
gamma_sdf   = funfitxy(fspace1, grid, sdf);
gamma_erk   = funfitxy(fspace1, grid, E_Rk);

muc_prime = E*gamma_lam;
sdf_prime = E*gamma_sdf;
ERk_prime =E*gamma_erk;

mul = 1 - beta*u.^gamma.*(gk).^(-gamma).*sdf_prime.*R./constr;
mu = max(0,mul);

pol_func.ERk = ERk_prime./Q;
pol_func.premium = ERk_prime./Q-R;
pol_func.mu = mu;
