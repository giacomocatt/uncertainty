%cd('C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\codes')
%addpath C:\dynare\5.3\matlab
var  
%%States
d z w theta nu

%%Controls
c i y l n u v ev psi %mu

%%Prices
Q R Rk exc_r_k
;

varexo eps_z eps_d eps_w;

parameters  

%%HH
beta gamma chi rho

%%Producers
delta alpha a1 csi zeta f

%%Banks
omega sigma theta_mean

%%Shocks
rho_z sd_z rho_w sd_w
 
%%Steady State
c_SS i_SS Q_SS R_SS y_SS d_SS l_SS n_SS;

alpha        = 0.33;
delta        = 0.025;
beta         = 0.995;
gamma        = 4;
rho          = 0.5;
csi          = 1;
z_mean       = 0;
l_SS          = 0.04;%%1/1.4;%%0.04;
y_SS        = l_SS^(1-alpha);
a1      = delta^(-csi);
theta_mean = 0.2191;
R_SS    = beta^(-1);
Q_SS    =1;
i_SS    = delta;
c_SS    = y_SS - i_SS;
zeta     =0.0933;
f = 0.4233;
Rk_SS   = (1-delta-zeta*f)*Q_SS + (alpha+f*(1-alpha))*y_SS;
chi = (1-f)*(1-alpha)*l_SS^(-alpha)+f*zeta*Q_SS/(l_SS);
omega = 0.002/4;
sigma = 0.972;%%1-0.24/12;
n_SS    = (sigma*(Rk_SS - R_SS*Q_SS)+omega*Q_SS)/(1-sigma*R_SS);
d_SS    = R_SS*(Q_SS-n_SS);
RN = (Rk_SS - d_SS)/n_SS;
psi_SS  =theta_mean*Q_SS/(beta*sigma*(Rk_SS-d_SS))-(1-sigma)/sigma;
mu_SS = 1-(1-sigma + sigma*psi_SS)*n_SS/(theta_mean*Q_SS);
rho_z        = 0.95;  
sd_z =0.01;
rho_w = 0.8;
sd_w=0.01;


model; 

%%HH
u = c-chi*l;

v = ((1-beta)*u^(1-rho) + beta*(ev)^((1-rho)))^(1/(1-rho));

ev^(1-gamma) = v(+1)^(1-gamma)*(1-delta+i);

beta * (u(+1)/u)^(-rho)*(ev/v(+1))^((gamma-rho))*(1-delta+i)^(-gamma)*R = 1;

chi*l= (1-f)*(1-alpha)*exp(z)*l^(1-alpha) + f*zeta*Q;

%%Banks

Rk = (1-delta-zeta*f)*Q + (alpha+f*(1-alpha))*y;

psi*n = theta*Q;

%mu = 1- exp(w(-1)*eps_d)*beta*(u(+1)/u)^(-rho)*(ev/v(+1))^((gamma-rho))*(1-delta+i)^(-gamma)*(1-sigma + sigma *psi(+1))*R*n/(theta*Q);

%beta*(u(+1)/u)^(-rho)*(ev/v(+1))^((gamma-rho))*(1-delta+i)^(-gamma)*(1-sigma + sigma *psi(+1))*(Rk(+1)/Q-R) = mu*theta;

psi*n = beta*(u(+1)/u)^(-rho)*(ev/v(+1))^((gamma-rho))*(1-delta+i)^(-gamma)*(1-sigma + sigma *psi(+1))*((Rk(+1)-d));

%%Producers

Q = a1*i^(csi);

y = exp(z)*l^(1-alpha);

%%Mkt clearing

c + i=y;%%a1*i^(1+csi)/(1+csi) - csi*delta/(1+csi) = y;

%%LOM

n= (sigma*(Rk - d(-1)) + omega*Q)/(1-delta+i);

d = R*(Q - n )*exp(nu);

z = rho_z*z(-1)+ eps_z;
log(theta) = log(theta_mean);
nu = (w(-1)*eps_d);
w = (1-rho_w)*sd_z+rho_w*w(-1) + eps_w;

exc_r_k =Rk/Q - R;


end;


initval;
theta = theta_mean;
%mu = mu_SS;
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
v = c_SS - chi*l_SS;
ev = (c_SS - chi*l_SS);
w = sd_z;
exc_r_k = Rk_SS/Q_SS - R_SS;
end;

steady(solve_algo=1, maxit=10000);
check;



shocks;
var eps_z;
stderr 0.01;
var eps_d;
stderr 0.1;
var eps_w;
stderr 1;
end;



stoch_simul(order=3 ,pruning, k_order_solver,irf=0);

IRF_periods=20;

burnin=5000; %periods for convergence

shock_mat_with_zeros=zeros(burnin+IRF_periods,M_.exo_nbr); %shocks set to 0 to simulate without uncertainty
IRF_no_shock_mat = simult_(M_,options_,oo_.dr.ys,oo_.dr,shock_mat_with_zeros,options_.order)'; %simulate series
stochastic_steady_state=IRF_no_shock_mat(1+burnin,:); % stochastic_steady_state/EMAS is any of the final points after burnin

shock_mat = zeros(burnin+IRF_periods,M_.exo_nbr);
shock_mat(1+burnin,strmatch('eps_w',M_.exo_names,'exact'))= 1;
IRF_mat = simult_(M_,options_,oo_.dr.ys,oo_.dr,shock_mat,options_.order)';

IRF_mat_percent_from_SSS = ((IRF_mat(1+burnin+1:1+burnin+IRF_periods,:))-(IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:)))./IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:); %only valid for variables not yet logged

y_pos 	= strmatch('y',M_.endo_names,'exact');
c_pos 	= strmatch('c',M_.endo_names,'exact');
inv_pos = strmatch('i',M_.endo_names,'exact');
l_pos 	= strmatch('l',M_.endo_names,'exact');
rk_pos 	= strmatch('Rk',M_.endo_names,'exact');
r_pos 	= strmatch('R',M_.endo_names,'exact');
erk_pos = strmatch('exc_r_k',M_.endo_names,'exact');
q_pos = strmatch('Q',M_.endo_names,'exact');
%%psi_pos    = strmatch('psi',M_.endo_names,'exact');
n_pos = strmatch('n',M_.endo_names,'exact');

y_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,y_pos);
c_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,c_pos);
inv_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,inv_pos);
n_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,n_pos);
rk_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,rk_pos);
r_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,r_pos);
erk_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,erk_pos);
q_vola_IRF= 100*IRF_mat_percent_from_SSS(:,q_pos);
%%psi_vola_IRF      = 100*IRF_mat_percent_from_SSS(:,psi_pos);
l_vola_IRF      = 100*IRF_mat_percent_from_SSS(:,l_pos);

hh=figure;
figure(hh)   
subplot(3,3,1)
hold on
plot(y_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$y$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,2)
hold on
plot(c_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$c$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)
%ylim([-0.3 0.1]);set(gca,'YTick',[-0.3:0.1:0.1],'FontSize',12);


figure(hh)   
subplot(3,3,3)
hold on
plot(inv_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$i$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)
%ylim([-0.6 0.4]);set(gca,'YTick',[-0.6:0.2:0.4],'FontSize',12);

figure(hh)   
subplot(3,3,4)
hold on
plot(l_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$\ell$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,5)
hold on
plot(erk_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Excess $R^k$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)


figure(hh)   
subplot(3,3,6)
hold on
plot(rk_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$R^k$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,7)
hold on
plot(r_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$R$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,8)
hold on
plot(q_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$Q$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,9)
hold on
plot(n_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('$n$','interpreter', 'latex','FontSize',14)
ylabel('Percent','FontSize',12)