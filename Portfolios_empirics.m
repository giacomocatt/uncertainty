cd 'C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty'
addpath("C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\Data")
addpath("C:\Users\giaco\Desktop\phd\Useful data")

close all
clear
controls_macro_data = readtable("Assets.xlsx", 'Sheet','controls_macro');
controls_fin_data = readtable("Assets.xlsx", 'Sheet','controls_fin');
asset_data = readtable("Assets.xlsx", 'Sheet','assets');
shock_data = readtable("Assets.xlsx", 'Sheet','shocks');
asset = table2array(asset_data(:, 3:6));
controls = [table2array(controls_macro_data(:, end-2:end)), table2array(controls_fin_data(:, [2:5]))];
capratio = table2array(controls_fin_data(:, end-2));
shock = table2array(shock_data(:, [2,4,5]));
for i =2:4
asset(:,i) = movmean(asset(:,i),12);
end

L = min([size(asset,1), size(controls,1), size(shock,1)]);
n = size(asset,2);
hor_back = 0;
hor =24;
p=[6,12,12,6];
q = [1,6,6,6];
shock_size = [std(shock(:,1)), std(shock(:,2)),std(shock(:,end))];
proj = zeros(3*n,hor+hor_back+1);
upproj95 = zeros(3*n,hor+hor_back+1);
lowproj95 = zeros(3*n,hor+hor_back+1);
upproj90 = zeros(3*n,hor+hor_back+1);
lowproj90 = zeros(3*n,hor+hor_back+1);

for i = 0:2
for j = 1:n
    lags = lagmatrix([asset(1:L,j)], q(j):p(j));
    xx = rmmissing([asset(p(j)+1:L,j), controls(p(j)+1:L,:), lags(p(j)+1:L,:),shock(p(j)+1:L,end-i)]);
    for h=1:hor
        %lags_shock = rmmissing(lagmatrix(xx(:,end-2:end-1),1:h));
        %zh =  [lags_shock(1:end,end-3), xx(1:end-h,2:end)];
        zh = [xx(1:end-h,2:end)];
        %zh = [xx(h:end-1,2:end-3), xx(1:end-h,end-2:end),xx(1:end-h,end).*xx(1:end-h,11)];
        x = xx(h+1:end,1);
        lp = fitlm(zh,x);
       proj(j+i*n,h+hor_back+1) = shock_size(end-i)*table2array(lp.Coefficients(end,1));
       upproj95(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) + shock_size(end-i)*table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
       lowproj95(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) -shock_size(end-i)* table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
       upproj90(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) + shock_size(end-i)*table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
       lowproj90(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) -shock_size(end-i)* table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));

        %proj(j+n,h+hor_back+1)      =               shock_size(2)*table2array(lp.Coefficients(end-2,1));
        %upproj95(j+n,h+hor_back+1) = proj(j+n,h+1) +shock_size(2)*table2array(lp.Coefficients(end-2,2))*tinv(0.975,size(zh,1)-size(zh,2));
        %lowproj95(j+n,h+hor_back+1) = proj(j+n,h+1)-shock_size(2)*table2array(lp.Coefficients(end-2,2))*tinv(0.975,size(zh,1)-size(zh,2));
        %upproj90(j+n,h+hor_back+1) = proj(j+n,h+1) +shock_size(2)*table2array(lp.Coefficients(end-2,2))*tinv(0.95,size(zh,1)-size(zh,2));
        %lowproj90(j+n,h+hor_back+1) = proj(j+n,h+1)-shock_size(2)*table2array(lp.Coefficients(end-2,2))*tinv(0.95,size(zh,1)-size(zh,2));
    end
end
end

newnames = ["SP500",'Bond Spread', "EBP", "GZ"];
%newnames = ["Portfolios", "SP500", "EBP", "GZ", "10yr-3m", "FX", "Interm"];
%newnames = ["Small-Low", "ME1-BM2", "Small-Hi", "Big-Low", "ME2-BM2", "Big-High"];
%newnames = ["Capital ratio ($\eta_t$)", "Risk factor ($(\eta_t - \rho_0 - \rho \eta_{t-1})/\eta_{t-1}$)"];
%newnames = ["Ret", "10yr Ret"];


color = [.7 .7 .7];


hor_back =1;
hor = hor-1;
fig1 = figure;
tcl = tiledlayout(3,4);
nexttile(tcl)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(1,:), fliplr(lowproj90(1,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(1,:), fliplr(lowproj95(1,:))], color, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(1,:), "MarkerFaceColor",color)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(1),'interpreter', 'latex')
ylabel('VXO', 'interpreter', 'latex')
for k=2:n
nexttile(tcl);
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k,:), fliplr(lowproj90(k,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k,:), fliplr(lowproj95(k,:))], color, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(k,:), "MarkerFaceColor",color)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(k),'interpreter', 'latex')
hold off
end
nexttile(tcl);
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(1+n,:), fliplr(lowproj90(1+n,:))],color, EdgeColor="none", FaceAlpha = 0.3)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(1+n,:), fliplr(lowproj95(1+n,:))], color, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(1+n,:), "MarkerFaceColor",color)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(1),'interpreter', 'latex')
ylabel('JLN financial', 'interpreter', 'latex')
for k=n+2:2*n
nexttile(tcl)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k,:), fliplr(lowproj90(k,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k,:), fliplr(lowproj95(k,:))], color, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(k,:), "MarkerFaceColor",color)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(k-n),'interpreter', 'latex')
hold off
end
nexttile(tcl);
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(1+2*n,:), fliplr(lowproj90(1+2*n,:))],color, EdgeColor="none", FaceAlpha = 0.3)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(1+2*n,:), fliplr(lowproj95(1+2*n,:))], color, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(1+2*n,:), "MarkerFaceColor",color)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(1),'interpreter', 'latex')
ylabel('JLN macro', 'interpreter', 'latex')
for k=2*n+2:3*n
nexttile(tcl)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k,:), fliplr(lowproj90(k,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k,:), fliplr(lowproj95(k,:))], color, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(k,:), "MarkerFaceColor",color)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(k-2*n),'interpreter', 'latex')
hold off
end
han=axes(fig1,'visible','off'); 
%han.Title.Visible='on';
%han.XLabel.Visible='on';
%han.YLabel.Visible='on';
%ylabel(han,'p.p.','interpreter', 'latex');
%xlabel(han,'Months','interpreter', 'latex');
print -depsc stocks_bonds_lin.eps