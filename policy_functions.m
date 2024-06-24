function [policy_functions] = policy_functions(coeffs, grid, fspace, M_, shocks, prob)

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

cons       = funeval(coeffs(:,1), fspace, grid);
R   = funeval(coeffs(:,2), fspace, grid);
%psi = funeval(coeffs(:,3), fspace, grid);

policy_functions.C = cons;

policy_functions.R = R;

%policy_functions.psi =psi;
%=========================================================================
%                COMPUTE ENDOGENOUS VARIABLES GIVEN GUESS
%=========================================================================

policy_functions.L     = ((1-alpha)/chi*exp(grid(:,2) + alpha*grid(:,1))./cons.^gamma).^(1/(alpha+nu));
     
policy_functions.Y       =  (policy_functions.L.^(1-alpha)).*exp(grid(:,2) + alpha*(grid(:,1)));

policy_functions.inve    =  max(policy_functions.Y-cons,0.01);

policy_functions.K_prime =  (1-delta)*exp(grid(:,1)) + policy_functions.inve;

policy_functions.Q       =  a1*(policy_functions.inve./exp(grid(:,1))).^(csi);

policy_functions.k_ret   =  (1-delta)*policy_functions.Q+(alpha*policy_functions.Y./exp(grid(:,1)));

policy_functions.N_net   = policy_functions.k_ret.*exp(grid(:,1)) - exp(grid(:,3));

policy_functions.N_prime = sigma*(policy_functions.N_net) + omega*policy_functions.Q.*policy_functions.K_prime;

policy_functions.D_prime = R.*(policy_functions.Q.*policy_functions.K_prime - policy_functions.N_prime);

%n=1;
policy_functions.constr = exp(theta_mean)*policy_functions.Q.*policy_functions.K_prime./policy_functions.N_prime;
%
policy_functions.psi     = ...%policy_functions.constr.*exp(n*policy_functions.constr)./(exp(n)+exp(n*policy_functions.constr)) + ones(N,1).*exp(n)./(exp(n)+exp(n*policy_functions.constr));
                            max(policy_functions.constr,1);

muc     = cons.^(-gamma);
sdf  = cons.^(-gamma).*(1-sigma+sigma*policy_functions.psi);
E_k_ret = sdf.*policy_functions.k_ret;

grid_prime   = [log(policy_functions.K_prime), rho_z*grid(:,2) + (1-rho_z)*z_mean, log(policy_functions.D_prime)];

gamma_lam   = funfitxy(fspace, grid, muc);
gamma_sdf   = funfitxy(fspace, grid, sdf);
gamma_erk   = funfitxy(fspace, grid, E_k_ret);

policy_functions.sdf_prime = expect_func(grid_prime, gamma_sdf, sd_z, shocks, prob, fspace);
Ek_ret_prime =expect_func(grid_prime, gamma_erk, sd_z, shocks, prob, fspace);

policy_functions.exc_ret_k = Ek_ret_prime./policy_functions.Q - R;
policy_functions.mu = max(1-beta*cons.^gamma.*policy_functions.sdf_prime.*R./policy_functions.constr,0);

end