
cd 'C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty'
addpath("C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\Data")
addpath("C:\Users\giaco\Desktop\phd\Useful data")

unc_shock_data = readtable("unc_shock_quarterly");
unc_shock = table2array(unc_shock_data(41:end,[2:4]));

%GROWTH RATES
clear all
close all
data = readtable("Annual_growth_real.xlsx");
X = table2array(data(:,2:end));
hor =10;
p=2;
proj = zeros(9,hor);
upproj95 = zeros(9,hor);
lowproj95 = zeros(9,hor);
upproj90 = zeros(9,hor);
lowproj90 = zeros(9,hor);

%Specification: x_it+h = const + u_t+h + ... + u_t + Sum_j x_jt-1
close all
Z = rmmissing(X);
Zlags = lagmatrix(Z(:,1:end-3), 1:p);
z = [Zlags(p+1:end,:), Z(p+1:end,end-2)];           %uf shock = Z(...,end), um shock = Z(...,end-2)
for j = 1:13
    for h=1:hor
        zh = z(1:end-h,:);
        zh(:,j) = [];
        x = Z(p+1+h:end,j);
        lp = fitlm(zh,x);
        proj(j,h) = table2array(lp.Coefficients(end,1));
        upproj95(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
    end
end
newnames= ["$\%\Delta$ Credit", "$\%\Delta$ Treasuries", "$\%\Delta$ Net Worth" , "$\%\Delta$ Private Safe Assets","$\%\Delta$ I", "$\%\Delta$ C",...
    "$\%\Delta$ Hours", "3M Real Rate", "10Y Real Rate", "$S\&P500$ $\log Q$", "$S\&P500$ $R$", "GZ EBP", "$S\&P500$ Excess $R$", ...
    "$\pi$", "$JLN^m$", "$JLN^f$", "Real GDP", "shock macro unc", "shock ip", "shock fin unc"];
fig = figure;
for k=1:4
    subplot(2,2,k)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,:), fliplr(lowproj95(k,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(k,:), fliplr(lowproj90(k,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(k,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc um_gr_finsector.eps

fig = figure;
for k=5:7
    subplot(1,3,k-4)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,:), fliplr(lowproj95(k,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(k,:), fliplr(lowproj90(k,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(k,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc uf_gr_realvars.eps

fig = figure;
for k=8:13
    subplot(2,3,k-7)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,:), fliplr(lowproj95(k,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(k,:), fliplr(lowproj90(k,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(k,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc um_gr_finvars.eps

%Bonds and Loans, Specification: x_it+h = const + u_t+h + ... + u_t + Sum_j x_jt-1
clear all
close all
data = readtable("Annual_growth_bonds_loans.xlsx");
X = table2array(data(:,2:end));
hor =10;
p=2;
proj = zeros(9,hor);
upproj95 = zeros(9,hor);
lowproj95 = zeros(9,hor);
upproj90 = zeros(9,hor);
lowproj90 = zeros(9,hor);

close all
Z = rmmissing(X);
Zlags = lagmatrix(Z(:,1:end-3), 1:p);
z = [Zlags(p+1:end,:), Z(p+1:end,end-2)];
for j = 1:12
    for h=1:hor
        zh = z(1:end-h,:);
        zh(:,j) = [];
        x = Z(p+1+h:end,j);
        lp = fitlm(zh,x);
        proj(j,h) = table2array(lp.Coefficients(end,1));
        upproj95(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
    end
end
newnames= ["$\%\Delta$ Corp. Bonds", "$\%\Delta$ Treasuries", "$\%\Delta$ Loans" ,"$\%\Delta$ I", "$\%\Delta$ C",...
    "$\%\Delta$ Hours", "3M Real Rate", "10Y Real Rate", "$S\&P500$ $\log Price$", "$S\&P500$ $R$", "GZ Bond Spread", "$S\&P500$ $R^{Exc}$", ...
    "$\pi$", "$JLN^m$", "$JLN^f$", "Real GDP", "shock macro unc", "shock ip", "shock fin unc"];
fig = figure;
subplot(1,2,1)
hold on
fill([1:10, fliplr(1:10)], [upproj95(1,:), fliplr(lowproj95(1,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(1,:), fliplr(lowproj90(1,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(1,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(1),'interpreter', 'latex')
xlim([1,10])
hold off
subplot(1,2,2)
hold on
fill([1:10, fliplr(1:10)], [upproj95(3,:), fliplr(lowproj95(3,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(3,:), fliplr(lowproj90(3,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(3,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(3),'interpreter', 'latex')
xlim([1,10])
hold off
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc um_gr_bl.eps

%Specification: x_it+h = const + u_t+h + ... + u_t + x_it-1 + x_it-1^2 + x_it-1^3 + gdp_t-1

close all
for j = 1:12
    Z = rmmissing(X(:,[j,end-3, end]));
    Zlags = lagmatrix(Z(:,1:end-1), 1:p);
    z = [Zlags(p+1:end,:), Zlags(p+1:end,:).^2, Zlags(p+1:end,:).^3, Z(p+1:end,end)];
    for h=1:hor
        zh = z(1:end-h,:);
        x = Z(p+1+h:end,1);
        lp = fitlm(zh,x);
        proj(j,h) = table2array(lp.Coefficients(end,1));
        upproj95(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
    end
end

newnames= ["$\Delta \log$ Credit", "$\Delta \log$ Treasuries", "$\Delta \log$ Equity", "$\Delta \log$ I", "$\Delta \log$ C", "$\Delta \log$ Hours", "3M Real Rate", ...
    "10Y Real Rate", "$\log Q^{S\&P500}$", "$R^{S\&P500}$", "GZ Bond Spread", "$S\&P500$ $R^{Exc}$", "$\pi$", "$JLN^m$", "$JLN^f$", "Real GDP", "shock macro unc", "shock ip", "shock fin unc"];
fig = figure;
for k=1:12
    subplot(4,3,k)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,:), fliplr(lowproj95(k,:))], [.9 .9 .9], EdgeColor="none")
hold on
plot(1:10, proj(k,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

%PERCENT DEVIATION FROM TREND
close all
clear all
data = readtable("Levels.xlsx");
X = table2array(data(:,2:end));
tobefiltered = X(:,[1,2,4,6,17]);
T = zeros(size(tobefiltered));
C = zeros(size(tobefiltered));
for i=1:5
    [T(:,i), C(:,i)] = hpfilter(tobefiltered(:,i));
end
[Tdep,Cdep] = hpfilter(X(:,3)); %GFCF
[Tinv,Cinv] = hpfilter(X(:,5)); %GFCF
X(:,[1,2,4,6,17]) = C(:,:)./T(:,:)*100;
X(9:end,3) = Cdep./Tdep*100;
X(5:end,5) = Cinv./Tinv*100;

controls = X(:,[9,15:17]);
hor =10;
p=2;
proj = zeros(9,hor);
upproj95 = zeros(9,hor);
lowproj95 = zeros(9,hor);
upproj90 = zeros(9,hor);
lowproj90 = zeros(9,hor);

Z = rmmissing(X);
Zlags = lagmatrix(Z(:,1:end-3), 1:p);
shock = 0;                                          %fin = 0, macro = 2
z = [Zlags(p+1:end,:), Z(p+1:end,end-shock)];
for j = 1:13
    for h=1:hor
        zh = z(1:end-h,:);
        zh(:,j) = [];
        x = Z(p+1+h:end,j);
        lp = fitlm(zh,x);
        proj(j,h) = table2array(lp.Coefficients(end,1));
        upproj95(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j,h) = proj(j,h) + table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j,h) = proj(j,h) - table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
    end
end

newnames= ["Credit", "Treasuries", "Fin. Sector Net Worth", "Pvt Safe Assets", "I", "C", "Hours", "3M Real Rate", ...
    "10Y Real Rate", "$S\&P500$ $\log Q$", "$S\&P500$ $R$", "GZ Bond Spread", "$S\&P500$ Exc. $R$", "$\pi$", "$JLN^m$", "$JLN^f$", ...
    "$\Delta \log$ Real GDP", "shock macro unc", "shock ip", "shock fin unc"];
fig = figure;
for k=1:4
    subplot(2,2,k)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,1:10), fliplr(lowproj95(k,1:10))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(k,1:10), fliplr(lowproj90(k,1:10))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(k,1:10), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc um_lv_finsector.eps

fig = figure;
for k=5:7
    subplot(1,3,k-4)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,:), fliplr(lowproj95(k,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(k,:), fliplr(lowproj90(k,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(k,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc uf_lv_realvars.eps

fig = figure;
for k=8:13
    subplot(2,3,k-7)
hold on
fill([1:10, fliplr(1:10)], [upproj95(k,:), fliplr(lowproj95(k,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(k,:), fliplr(lowproj90(k,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(k,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title(newnames(k),'interpreter', 'latex')
xlim([1,10])
hold off
end
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc um_lv_finvars.eps

fig = figure;
    subplot(1,2,1)
hold on
fill([1:10, fliplr(1:10)], [upproj95(1,:), fliplr(lowproj95(1,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(1,:), fliplr(lowproj90(1,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(1,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title('Corp. Bonds','interpreter', 'latex')
xlim([1,10])
hold off

subplot(1,2,2)
hold on
fill([1:10, fliplr(1:10)], [upproj95(3,:), fliplr(lowproj95(3,:))], [.9 .9 .9], EdgeColor="none")
hold on
fill([1:10, fliplr(1:10)], [upproj90(3,:), fliplr(lowproj90(3,:))], [.8 .8 .8], EdgeColor="none")
hold on
plot(1:10, proj(3,:), "k")
hold on
plot(zeros(1,10), 'r')
box on
grid on
title('Loans','interpreter', 'latex')
xlim([1,10])
hold off

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Quarters','interpreter', 'latex');

print -depsc uf_lv_bl.eps