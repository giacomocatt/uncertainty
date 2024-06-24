%addpath C:\dynare\5.3\matlab
var  
%%States
d z %%w

%%Controls
c i y l n psi u

%%Prices
Q R Rk
;

varexo eps_z eps_d;%%eps_w;

parameters  

%%HH
beta gamma chi

%%Producers
delta alpha a1 csi zeta f

%%Banks
omega sigma theta_mean

%%Shocks
rho_z sd_z rho_w sd_w
 
%%Steady State
c_SS i_SS Q_SS R_SS y_SS d_SS l_SS n_SS;

alpha        = 0.33;
delta        = 0.05;
beta         = 0.995;
gamma        = 1;
csi          = 1;
z_mean       = 0;
l_SS          = 0.1;%%1/1.4;%%0.1;
y_SS        = l_SS^(1-alpha);
a1      = delta^(-csi);
theta_mean = 0.31;
R_SS    = beta^(-1);
Q_SS    = 1;
chi= (1-f)*(1-alpha)*l_SS^(-alpha) + f*zeta*Q_SS/l_SS;
%%y_SS    = (spread-(1-delta) + R_SS)/alpha;
i_SS    = delta;
c_SS    = y_SS - i_SS;
Rk_SS   = (1-delta)*Q_SS+alpha*y_SS;
zeta     =0.0933;
f = 0.4233;
chi = (1-f)*(1-alpha)*l_SS^(-alpha)+f*zeta*Q_SS/(l_SS);
omega = 0.002/12;
sigma = 0.97;%%1-0.24/12;
n_SS    = (sigma*(Rk_SS/Q_SS - R_SS)+omega)*Q_SS/(1-sigma*R_SS);
d_SS    = R_SS*(Q_SS-n_SS);
RN = (Rk_SS - d_SS)/n_SS;
psi_SS  =(theta_mean*Q_SS/(beta*(Rk_SS-d_SS))-(1-sigma))/sigma;
rho_z        = 0.9;  
sd_z =0.1;
rho_w = 0.8;
sd_w=0.1;


model; 

%%HH
u = c-chi*l;

beta * (u(+1)/u)^(-gamma)*(1-delta+i)^(-gamma)*R = 1;

chi*l= (1-f)*(1-alpha)*exp(z)*l^(1-alpha) + f*zeta*Q;

%%Banks

Rk = (1-delta)*Q + alpha*y;

psi*n= theta_mean*Q;

%%Producers

Q = a1*i^(csi);

y = exp(z)*l^(1-alpha);

%%Mkt clearing

c + a1*i^(1+csi)/(1+csi) - csi*delta/(1+csi) = y;

%%LOM

n= (sigma*(Rk - d(-1)*exp(eps_d)) + omega*Q)/(1-delta+i);

d = R*(Q - n );

psi = beta * (u(+1)/u)^(-gamma)*(1-delta+i)^(-gamma)*(1-sigma + sigma *psi(+1))*(Rk(+1) - d)/n;

z = rho_z*z(-1)+ eps_z;
%%w = (1-rho_w)*sd_z+rho_w*w(-1) + eps_w;


end;


initval;
R   = R_SS;  
Q   = Q_SS;    
i   = i_SS;   
c   = c_SS ;
l=l_SS;
psi = psi_SS;
d = d_SS;
n = n_SS;
u = c_SS - chi*l_SS;
y =y_SS;
l=l_SS;
Rk = Rk_SS;
end;

steady;
check;



shocks;
var eps_z;
stderr sd_z;
var eps_d;
stderr 0.1;
end;



stoch_simul(order=2 , irf=0);