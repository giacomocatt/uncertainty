function [R_and_psi] = R_psi_dz(policy, policy_next,grid,fspace, M_, shocks, prob)
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

R_t_1 = policy_next(N+1:2*N);
psi_t_1 = policy(2*N+1:3*N);
l        =  exp(policy_next(1:N));

Q       = (chi*l - (1-f)*(1-alpha)*exp(grid(:,2)).*l.^(1-alpha))/(zeta*f);

y       = exp(grid(:,2)).*(l).^(1-alpha);
       
i       = (Q/a1).^(1/csi);

c       =  max(y-i,0.01);

u       = c - chi*l;

Rk      =  (1-delta)*Q+alpha*y;

N_net  = max(Rk - exp(grid(:,1)),0.1);

N_prime= (sigma*(N_net) + omega*Q)./(1-delta+i);

D_prime = R_t_1.*(Q - N_prime);

muc = u.^(-gamma);

gamma_muc = funfitxy(fspace, grid, muc);
grid_prime = [log(D_prime),  rho_z*grid(:,2) + (1-rho_z)*z_mean];
E = expect_coeff_dz(grid_prime,unc, sd_z, shocks, prob, fspace);
Emuc = E*gamma_muc;
sdf = muc.*(psi_t_1);
gamma_sdf = funfitxy(fspace, grid, sdf);
Esdf = E*gamma_sdf;
%disc_Rk = muc.*(1-sigma + sigma*psi_t_1).*Rk;
%gamma_rk = funfitxy(fspace, grid, disc_Rk);
%ERk = E*gamma_rk;

R = max(min((1-delta+i).^gamma./(beta*u.^gamma.*Emuc),2),0.1);

unconstr = ...%beta*u.^gamma.*(1-delta+i).^(-gamma).*R.*Esdf;%...
           1-sigma+sigma* Esdf./Emuc;
constr = theta_mean*Q./(N_prime);
psi     = max(constr, unconstr);
 %beta*u.^gamma.*(1-delta+i).^(-gamma).*(ERk -  ((1-sigma)*Emuc+sigma* Esdf).*D_prime)./N_net;
R_and_psi = [R;psi];

end