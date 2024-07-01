function [pol_func] = policy_functions_dz(coeffs, GG, params, fspace, inve)
beta        = params(1);       
gamma       = params(2);  
chi         = params(3);
delta       = params(4);  
alpha       = params(5);       
a1          = params(6);                  
csi         = params(7);   
zeta        = params(8); 
f           = params(9);
omega       = params(10);
sigma       = params(11);
z_mean      = params(12);
sd_z        = params(13); 
rho_z       = params(16);
theta_mean  = params(end-1);
unc         = params(end);
N          = size(GG,1);
policy(1:N) = funeval(coeffs(:,1), fspace, GG);
policy(N+1:2*N)   = funeval(coeffs(:,2), fspace, GG);
policy(2*N+1:3*N)   = funeval(coeffs(:,3), fspace, GG);

if size(policy,1) ==1
    l        =  exp(policy(1:N))';
    R           =  policy(N+1:2*N)';
    psi     = policy(2*N+1:3*N)';
else
    l        =  exp(policy(1:N));
    R           =  policy(N+1:2*N);
    psi     = policy(2*N+1:3*N);
end



Q = (chi*l - (1-f)*(1-alpha)*exp(GG(:,2)).*l.^(1-alpha))/(zeta*f);
y       = exp(GG(:,2)).*(l).^(1-alpha);
       
i       = (Q/a1).^(1/csi);

if inve == 0
    gk = 1-delta + i;
    cost =  a1*i.^(1+csi)/(1+csi)-csi*delta/(1+csi);
else
    gk = 1-delta + a1^(-1)*i.^(1-csi)/(1-csi)-csi*delta/(1-csi);
    cost = i;
end

c       =  max(y-cost,1e-6);

Rk      =  (1-delta-f*zeta)*Q+alpha*y +f*(1-alpha)*y;

n_net  = max(Rk - exp(GG(:,1)),0);

n_prime = (sigma*(n_net) + omega*Q)./(gk);

d_prime = max(R.*(Q - n_prime),1e-6);

u       = c - chi*l;

constr = theta_mean*Q./(n_prime);

pol_func.l = l;
pol_func.R = R;
pol_func.psi = psi;
pol_func.Q = Q;
pol_func.y = y;
pol_func.c = c;
pol_func.i = i;
pol_func.rk = Rk./Q;
pol_func.n_net = n_net;
pol_func.n_prime = n_prime;
pol_func.d_prime = d_prime;
pol_func.u = u;

muc     = u.^(-gamma);
sdf     = (1-sigma+sigma*psi).*muc;
E_Rk    = Rk;

pol_func.GG_prime   = [log(d_prime), rho_z*GG(:,2)+(1-rho_z)*z_mean];
%E = expect_coeff_dz(GG_prime,unc, sd_z, shocks, prob, fspace1);
%pol_func.E = E;
%
%gamma_lam   = funfitxy(fspace1, GG, muc);
%gamma_sdf   = funfitxy(fspace1, GG, sdf);
%gamma_erk   = funfitxy(fspace1, GG, E_Rk);
%
%muc_prime = E*gamma_lam;
%sdf_prime = E*gamma_sdf;
%ERk_prime =E*gamma_erk;
%
%mul = 1 - beta*u.^gamma.*(gk).^(-gamma).*sdf_prime.*R./constr;
%mu = max(0,mul);
%
%pol_func.ERk = ERk_prime./Q;
%pol_func.premium = ERk_prime./Q-R;
%pol_func.mul = mul;
%pol_func.mu = mu;
