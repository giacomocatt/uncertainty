function [residuals] = EE_dz(policy_next, policy, grid, fspace, M_, shocks, prob, inve)

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
residuals = zeros(N,1);

l        =  exp(policy_next(1:N));
%psi      = policy_next(N+1:2*N);
%=========================================================================
%                COMPUTE ENDOGENOUS VARIABLES GIVEN GUESS
%=========================================================================
Q       = (chi*l - (1-f)*(1-alpha)*exp(grid(:,2)).*l.^(1-alpha))/(zeta*f);

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

u       = c - chi*l;

Rk      =  gk.*Q+alpha*y -cost ;

n_net  = max(Rk - exp(grid(:,1)),0.1);

n_prime = (sigma*(n_net) + omega*Q)./(gk);

constr = theta_mean*Q./(n_prime);
%=========================================================================
%                         GUESSES FOR FUTURE VALUES
%=========================================================================

l_t_1        =  exp(policy(1:N));
R_t_1           =  policy(N+1:2*N);
psi_t_1     = policy(2*N+1:3*N);

Q_t_1 = (chi*l_t_1 - (1-f)*(1-alpha)*exp(grid(:,2)).*l_t_1.^(1-alpha))/(zeta*f);
y_t_1       = exp(grid(:,2)).*(l_t_1).^(1-alpha);
       
i_t_1       = (Q_t_1/a1).^(1/csi);

if inve == 0
    gk_t_1 = 1-delta + i_t_1;
    cost_t_1 =  a1*i_t_1.^(1+csi)/(1+csi)-csi*delta/(1+csi);
else
    gk_t_1 = 1-delta + a1^(-1)*i_t_1.^(1-csi)/(1-csi)-csi*delta/(1-csi);
    cost_t_1 = i_t_1;
end

c_t_1       =  max(y_t_1-cost_t_1,0.01);

Rk_t_1      =  gk_t_1.*Q_t_1+alpha*y_t_1 -cost_t_1;

n_net_t_1  = max(Rk_t_1 - exp(grid(:,1)),0.1);

n_prime_t_1 = (sigma*(n_net_t_1) + omega*Q_t_1)./(gk_t_1);

d_prime = R_t_1.*(Q_t_1 - n_prime_t_1);

u_t_1       = c_t_1 - chi*l_t_1;

%=========================================================================
%                         COMPUTE EXPECTATIONS
%=========================================================================

muc     = u_t_1.^(-gamma);
sdf     = (1-sigma+sigma*psi_t_1).*muc;
E_Rk    = sdf.*Rk_t_1;

grid_prime   = [log(d_prime), rho_z*grid(:,2)+(1-rho_z)*z_mean];
E = expect_coeff_dz(grid_prime,unc, sd_z, shocks, prob, fspace);

gamma_sdf   = funfitxy(fspace, grid, sdf);
gamma_erk   = funfitxy(fspace, grid, E_Rk);

sdf_prime = E*gamma_sdf;
ERgk =E*gamma_erk;

%=========================================================================
%                         EULER EQUATION
%=========================================================================
mul = 1 - beta*u.^gamma.*(gk).^(-gamma).*sdf_prime.*R_t_1./constr;
%mu = mul.*exp(x*mul)./(1+exp(x*mul));
mu = max(0,mul).^2;

RHS = beta*u.^gamma.*(gk).^(-gamma).*(ERgk./Q - sdf_prime.*R_t_1);

residuals(:,1) = RHS.^2 - theta_mean^2*mu;

%residuals(:,2) = 1-sigma+...
%    beta*sigma*beta*u.^gamma.*(1-delta+i).^(-gamma).*(ERgk - sdf_prime.*exp(d_prime))./n_prime - psi;

residuals        = reshape(residuals,[],1);

end