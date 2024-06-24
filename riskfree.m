function [coeff_r] = riskfree(coeffs, grid, fspace, M_, shocks, prob)

beta        = M_.params(1);       
gamma       = M_.params(2);  
rho_z       = M_.params(7);           
z_mean      = M_.params(8); 
sd_z        = M_.params(14); 
N = size(grid,1);
C     = coeffs;
MUC = C.^(-gamma);
grid_prime   = [grid(:,1),rho_z*grid(:,2) + (1-rho_z)*z_mean, grid(:,3)];
E = expect_coeff_tfp(grid_prime, sd_z, shocks, prob, fspace);
 
gamma_muc  = funfitxy(fspace, grid, MUC);

EMUC = E*gamma_muc;

coeff_r = (beta*C.^gamma.*EMUC).^(-1);
%coeff_r = funfitxy(fspace, grid, R);

end