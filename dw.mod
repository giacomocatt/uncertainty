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
omega sigma

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
y_w          =2.4;
l_SS         = y_w^(1-alpha);
a1      = delta^(-csi);
R_SS    = beta^(-1);
y_SS    = l_SS^(1-alpha);
Q_SS    = 1;
i_SS    = delta;
c_SS    = y_SS - i_SS;
Rk_SS   = R_SS;
zeta     =0.1;
f = 0.5;
chi = (1-f)*(1-alpha)*l_SS^(-alpha)+f*zeta*Q_SS/(l_SS);
omega = 0.002/12;
sigma =1-0.24/12;
n_SS    = omega*Q_SS/(1-sigma*R_SS);
d_SS    = R_SS*(Q_SS-n_SS);
psi_SS  =beta*(1-sigma)/(1-beta*sigma);
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

beta * (u(+1)/u)^(-gamma)*(1-delta+i)^(-gamma)* exp(-gamma*z(+1))*Rk(+1)/Q =1;

%%Producers

Q = a1*i^(csi);

y = l^(1-alpha);

%%Mkt clearing

c + i = y;

%%LOM

n= sigma*(Rk - d(-1)) + omega*Q;

d = R*(Q - n/(1-delta+i)*exp(-z(+1)) );

psi = beta * (u(+1)/u)^(-gamma)*(1-delta+i)^(-gamma)*exp(-gamma*z(+1))* (1-sigma + sigma *psi(+1))*R;

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
end;

steady;
check;



shocks;
var eps_z;
stderr sd_z;
var eps_w;
stderr sd_w;
end;



stoch_simul(order=3, pruning, irf=20);
IRF_periods=20;

burnin=5000; %periods for convergence

shock_mat_with_zeros=zeros(burnin+IRF_periods,M_.exo_nbr); %shocks set to 0 to simulate without uncertainty
IRF_no_shock_mat = simult_(M_,options_,oo_.dr.ys,oo_.dr,shock_mat_with_zeros,options_.order)'; %simulate series
stochastic_steady_state=IRF_no_shock_mat(1+burnin,:); % stochastic_steady_state/EMAS is any of the final points after burnin

shock_mat = zeros(burnin+IRF_periods,M_.exo_nbr);
shock_mat(1+burnin,strmatch('eps_w',M_.exo_names,'exact'))= 0.1;
IRF_mat = simult_(M_,options_,oo_.dr.ys,oo_.dr,shock_mat,options_.order)';

IRF_mat_percent_from_SSS = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,:)-IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:))./repmat(stochastic_steady_state,IRF_periods,1); %only valid for variables not yet logged

y_pos 	= strmatch('y',M_.endo_names,'exact');
c_pos 	= strmatch('c',M_.endo_names,'exact');
inv_pos = strmatch('i',M_.endo_names,'exact');
n_pos 	= strmatch('l',M_.endo_names,'exact');
rk_pos 	= strmatch('Rk',M_.endo_names,'exact');
r_pos 	= strmatch('R',M_.endo_names,'exact');
q_pos = strmatch('Q',M_.endo_names,'exact');
psi_pos    = strmatch('psi',M_.endo_names,'exact');

y_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,y_pos);
c_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,c_pos);
inv_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,inv_pos);
n_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,n_pos);
rk_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,rk_pos);
r_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,r_pos);
q_vola_IRF= 100*IRF_mat_percent_from_SSS(:,q_pos);
psi_vola_IRF      = 100*IRF_mat_percent_from_SSS(:,psi_pos);

hh=figure;
figure(hh)   
subplot(3,3,1)
hold on
plot(y_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Output','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,2)
hold on
plot(c_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Consumption','FontSize',14)
ylabel('Percent','FontSize',12)
%ylim([-0.3 0.1]);set(gca,'YTick',[-0.3:0.1:0.1],'FontSize',12);


figure(hh)   
subplot(3,3,3)
hold on
plot(inv_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Investment','FontSize',14)
ylabel('Percent','FontSize',12)
%ylim([-0.6 0.4]);set(gca,'YTick',[-0.6:0.2:0.4],'FontSize',12);

figure(hh)   
subplot(3,3,4)
hold on
plot(n_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Labor','FontSize',14)
ylabel('Percent','FontSize',12)


figure(hh)   
subplot(3,3,5)
hold on
plot(rk_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Rk','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,6)
hold on
plot(r_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('R','FontSize',14)
ylabel('Percent','FontSize',12)

figure(hh)   
subplot(3,3,7)
hold on
plot(q_vola_IRF,'b-','LineWidth',3)
plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);
title('Q','FontSize',14)
ylabel('Percent','FontSize',12)
