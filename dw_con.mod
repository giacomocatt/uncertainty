%addpath C:\dynare\5.3\matlab
var  
%%States
d w

%%Controls
c i y l n psi u z

%%Prices
Q R Rk
;

varexo eps_z eps_w;

parameters  

%%HH
beta gamma chi

%%Producers
delta alpha a1 csi zeta f

%%Banks
omega sigma theta_mean

%%Shocks
sd_z sd_w rho_w
 
%%Steady State
c_SS i_SS Q_SS R_SS y_SS d_SS l_SS n_SS;

alpha        = 0.33;
delta        = 0.05;
beta         = 0.997;
gamma        = 1;
csi          = 1;
z_mean       = 0;
y_SS          = 1;
spread = 0.01;
a1      = delta^(-csi);
R_SS    = beta^(-1);
Q_SS    = (alpha*y_SS)/(R_SS+spread - 1+delta);
%%y_SS    = (spread-(1-delta) + R_SS)/alpha;
l_SS         = y_SS^(1/(1-alpha));
i_SS    = delta;
c_SS    = y_SS - i_SS;
Rk_SS   = (1-delta)*Q_SS+alpha*y_SS;
zeta     =0.3;
f = 0.5;
chi = (1-f)*(1-alpha)*l_SS^(-alpha)+f*zeta*Q_SS/(l_SS);
omega = 0.002/12;
sigma = 0.90;%%1-0.24/12;
n_SS    = (sigma*spread+omega)*Q_SS/(1-sigma*R_SS);
d_SS    = R_SS*(Q_SS-n_SS);
psi_SS  =beta*(1-sigma)*(Rk_SS-d_SS)/n_SS/(1-beta*sigma*(Rk_SS-d_SS)/n_SS);
theta_mean = psi_SS*n_SS/Q_SS;
rho_w        = 0.8;  
sd_z =0.1;
sd_w=0.1;


model; 

%%HH
u = c-chi*l;

beta * (u(+1)/u)^(-gamma)*(1-delta+i)^(-gamma)* exp(-gamma*z(+1))*R = 1;

chi= (1-f)*(1-alpha)*l^(-alpha) + f*zeta*Q/l;

%%Banks

Rk = (1-delta)*Q + alpha*y;

psi*n/(1-delta+i)= theta_mean*Q;

%%Producers

Q = a1*i^(csi);

y = l^(1-alpha);

%%Mkt clearing

c + i = y;

%%LOM

n= sigma*(Rk - d(-1)) + omega*Q;

d = R*(Q - n/(1-delta+i)*exp(-z(+1)) );

psi = beta * (u(+1)/u)^(-gamma)*(1-delta+i)^(-gamma)*exp(-gamma*z(+1))*(1-sigma + sigma *psi(+1))*(Rk(+1) - d)/n;

z = exp(w(-1))*eps_z;
w = (1-rho_w)*sd_z + rho_w*w(-1) + eps_w;


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
w = sd_z;
end;

steady;
check;



shocks;
var eps_z;
stderr sd_z;
var eps_w;
stderr sd_w;
end;



stoch_simul(order=3, pruning, irf=0);