function [residuals] = EE_banks(policy_next, policy, grid, fspace, M_, shocks, prob)

beta        = M_.params(1);       
gamma       = M_.params(2);  
chi         = M_.params(4);  
nu          = M_.params(3);
delta       = M_.params(5);  
alpha       = M_.params(6);       
rho_z       = M_.params(7);           
z_mean      = M_.params(8);       
a1          = M_.params(9);                  
csi         = M_.params(10);   
vartheta    = M_.params(11); 
theta_mean  = M_.params(end);
omega       = M_.params(12);
sigma       = M_.params(13);
sd_z        = M_.params(14); 
N          = size(grid,1);
residuals = zeros(N,1);

%policy_next(1:N)       = funeval(coeff_next(1:N), fspace, grid)';
%policy_next(N+1:2*N)   = funeval(coeff_next(N+1:2*N), fspace, grid)';
%policy_next(2*N+1:3*N) = funeval(coeff_next(2*N+1:3*N), fspace, grid)';
%policy_next(3*N+1:4*N) = funeval(coeff_next(3*N+1:4*N), fspace, grid)';

cons        =  policy_next(1:N);
%R           =  policy_next(N+1:2*N);
%R = beta^(-1);

%=========================================================================
%                COMPUTE ENDOGENOUS VARIABLES GIVEN GUESS
%=========================================================================

L     = ((1-alpha)/chi*exp(grid(:,2) + alpha*grid(:,1))./cons.^gamma).^(1/(alpha+nu));
       
Y       =  (L.^(1-alpha)).*exp(grid(:,2) + alpha*(grid(:,1)));

inve    =  max(Y-cons,0.01);

K_prime =  (1-delta)*exp(grid(:,1)) + inve;

Q       =  a1*(inve./exp(grid(:,1))).^(csi);

k_ret   =  (1-delta)*Q+(alpha*Y./exp(grid(:,1)));

N_net   = max(k_ret.*exp(grid(:,1)) - exp(grid(:,3)),0.1);

N_prime = sigma*(N_net) + omega*Q.*K_prime;

constr     = exp(theta_mean)*Q.*K_prime./N_prime;

%=========================================================================
%                         GUESSES FOR FUTURE VALUES
%=========================================================================
%policy(1:N)       = funeval(coeffs(1:N), fspace, grid)';
%policy(N+1:2*N)   = funeval(coeffs(N+1:2*N), fspace, grid)';
%policy(2*N+1:3*N) = funeval(coeffs(2*N+1:3*N), fspace, grid)';

cons_t_1        =  policy(1:N);
R_t_1           =  policy(N+1:2*N);
%R_t_1           =  beta^(-1);

L_t_1     = ((1-alpha)/chi*exp(grid(:,2) + alpha*grid(:,1))./cons_t_1.^gamma).^(1/(alpha+nu));

Y_t_1       =  (L_t_1.^(1-alpha)).*exp(grid(:,2) + alpha*(grid(:,1)));

inve_t_1    =  max(Y_t_1-cons_t_1,0.01);

Q_t_1       =  a1*(inve_t_1./exp(grid(:,1))).^(csi);

K_prime_t_1 =  (1-delta)*exp(grid(:,1)) + inve_t_1;

k_ret_t_1   =  (1-delta)*Q_t_1+(alpha*Y_t_1./exp(grid(:,1)));

N_net_t_1   = max(k_ret_t_1.*exp(grid(:,1)) - exp(grid(:,3)),0.1);

N_prime_t_1 = sigma*(N_net_t_1) + omega*Q_t_1.*K_prime_t_1;

D_prime_t_1 = R_t_1.*(Q_t_1.*K_prime_t_1 - N_prime_t_1);

%constr = exp(theta_mean)*Q_t_1.*K_prime_t_1./N_prime_t_1;
%n = 5;
%psi_t_1     = constr.*exp(n*constr)./(exp(n*policy(2*N+1:3*N))+exp(n*constr))...
%            + policy(2*N+1:3*N).*exp(n*policy(2*N+1:3*N))./(exp(n*policy(2*N+1:3*N))+exp(n*constr));
psi_t_1     = policy(2*N+1:3*N);

%=========================================================================
%                         COMPUTE EXPECTATIONS
%=========================================================================

muc     = cons_t_1.^(-gamma);
sdf     = muc.*(1-sigma+sigma*psi_t_1);
E_k_ret = sdf.*k_ret_t_1;

grid_prime   = [log(K_prime_t_1),  rho_z*grid(:,2) + (1-rho_z)*z_mean, log(D_prime_t_1)];
E = expect_coeff_tfp(grid_prime, sd_z, shocks, prob, fspace);

gamma_lam   = funfitxy(fspace, grid, muc);
gamma_sdf   = funfitxy(fspace, grid, sdf);
gamma_erk   = funfitxy(fspace, grid, E_k_ret);

muc_prime = E*gamma_lam;
sdf_prime = E*gamma_sdf;
Ek_ret_prime =E*gamma_erk;

%muc_prime = expect_func(grid_prime, gamma_lam, sd_z, shocks, prob, fspace);
%sdf_prime = expect_func(grid_prime, gamma_sdf, sd_z, shocks, prob, fspace);
%Ek_ret_prime =expect_func(grid_prime, gamma_erk, sd_z, shocks, prob, fspace);

%=========================================================================
%                         COMPUTE IMPLIED PF
%=========================================================================
n = 1;
mul = 1 - beta*cons.^gamma.*sdf_prime.*R_t_1./constr;
mu = mul.*exp(n*mul)./(1+exp(n*mul));
%mu = max(mul,0);

residuals(:,1) = beta*cons.^gamma.*(Ek_ret_prime./Q - sdf_prime.*R_t_1) - exp(theta_mean)*mu;

%residuals(:,2) = chi*L.^(1+nu).*cons.^gamma - vartheta*Q.*exp(grid(:,1));

norm(residuals)^2;

residuals        = reshape(residuals,[],1);

end