function [coeff_r] = riskfree_dz(policy, grid, fspace, M_, shocks, prob, inve)

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
N = size(grid,1);
l = exp(policy);
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

C     = max(y-cost,0.01);
u = C - chi*l;
MUC = u.^(-gamma);
grid_prime   = [grid(:,1), rho_z*grid(:,2) + (1-rho_z)*z_mean];
E = expect_coeff_tfp(grid_prime, sd_z, shocks, prob, fspace);
 
gamma_muc  = funfitxy(fspace, grid, MUC);

EMUC = E*gamma_muc;

coeff_r = (gk).^gamma./(beta*u.^gamma.*EMUC);
%coeff_r = funfitxy(fspace, grid, R);

end