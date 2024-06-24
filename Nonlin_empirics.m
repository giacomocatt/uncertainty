cd 'C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty'
addpath("C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\Data")
addpath("C:\Users\giaco\Desktop\phd\Useful data")

close all
clear
controls_macro_data = readtable("Assets.xlsx", 'Sheet','controls_macro');
controls_fin_data = readtable("Assets.xlsx", 'Sheet','controls_fin');
asset_data = readtable("Assets.xlsx", 'Sheet','assets');
shock_data = readtable("Assets.xlsx", 'Sheet','shocks');
asset = table2array(asset_data(:, [3,5,6]));
%controls = [table2array(controls_macro_data(:, end-2:end)) table2array(controls_fin_data(:, [2:5]))];%vola
controls = [table2array(controls_macro_data(:, [end-6:end-3,end])), table2array(controls_fin_data(:, [2:5]))]; %JLN
capratio = table2array(controls_fin_data(:, end-1:end));
shock = table2array(shock_data(:, [2,4,5]));
for i =2:3
asset(:,i) = movmean(asset(:,i),12);
end

L = min([size(asset,1), size(controls,1), size(shock,1), size(capratio,1)]);
n = size(asset,2);
hor_back = 0;
hor =24;
%p = [12,6,6,6];
%q = [12,6,6,6];
p=[6,12,12,6];%vola
q = [1,6,6,6];%vola
select_control = [1,2,2,2];
shock_size = [std(shock(:,1)), std(shock(:,2)), std(shock(:,3))];
proj = zeros(2*n,hor+hor_back+1);
upproj95 = zeros(2*n,hor+hor_back+1);
lowproj95 = zeros(2*n,hor+hor_back+1);
upproj90 = zeros(2*n,hor+hor_back+1);
lowproj90 = zeros(2*n,hor+hor_back+1);
cr1 = capratio(:,1);
cr2 = capratio(:,2);
r_mean = [mean(cr1(cr1>=0)), mean(cr2(cr2>=0))];%mean(rmmissing(controls(:,end)));
rr = [mean(cr1(cr1<0)), mean(cr2(cr2<0))];%mean(rmmissing(controls(:,end))) - std(rmmissing(controls(:,end)));
%capratio(:,1) = double(capratio(:,1)<0);
capratio = lagmatrix(capratio,1);%lagmatrix(double(capratio(:,2)<0),2);

for j = 1:n
    lags = lagmatrix([asset(1:L,j), controls(1:L,:)], q(j):p(j));
    xx = rmmissing([asset(p(j)+1:L,j), controls(p(j)+1:L,:),lags(p(j)+1:L,:),...
        capratio(p(j)+1:L,select_control(j)), shock(p(j)+1:L,1), ...
        capratio(p(j)+1:L,select_control(j)).*shock(p(j)+1:L,1)]);
    for h=1:hor
        %lags_shock = rmmissing(lagmatrix(xx(:,end-2:end-1),1:h));
        %zh =  [lags_shock(1:end,end-3), xx(1:end-h,2:end)];
        zh = [xx(1:end-h,2:end)];
        %zh = [xx(h:end-1,2:end-3), xx(1:end-h,end-2:end),xx(1:end-h,end).*xx(1:end-h,11)];
        x = xx(h+1:end,1);
        lp = fitlm(zh,x);
       aa = shock_size(1)*table2array(lp.Coefficients(end,1));
       
        bb      =               shock_size(1)*table2array(lp.Coefficients(end-1,1));
       
        proj(j,h+hor_back+1)      =     bb+aa*r_mean(select_control(j));
        varproj = r_mean(select_control(j))^2*table2array(lp.Coefficients(end,2))^2 + table2array(lp.Coefficients(end-1,2))^2 ...
                +2*lp.CoefficientCovariance(end,end-1)*r_mean(select_control(j));
        upproj95(j,h+hor_back+1) = proj(j,h+1) ...
            + shock_size(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j,h+hor_back+1) = proj(j,h+1) ...
            - shock_size(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j,h+hor_back+1) = proj(j,h+1)...
            + shock_size(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j,h+hor_back+1) = proj(j,h+1)...
            - shock_size(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));

        proj(j+n,h+hor_back+1)      =     bb+aa*rr(select_control(j));
        varproj = rr(select_control(j))^2*table2array(lp.Coefficients(end,2))^2 + table2array(lp.Coefficients(end-1,2))^2 ...
                +2*lp.CoefficientCovariance(end,end-1)*rr(select_control(j));
        upproj95(j+n,h+hor_back+1) = proj(j+n,h+1) ...
            + shock_size(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j+n,h+hor_back+1) = proj(j+n,h+1) ...
            - shock_size(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j+n,h+hor_back+1) = proj(j+n,h+1)...
            + shock_size(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j+n,h+hor_back+1) = proj(j+n,h+1)...
            - shock_size(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));
    end
end

color1 =  [0.4940 0.1840 0.5560];%[0 0.4470 0.7410];%[0,0,1];
color2 = [0.3010 0.7450 0.9330];%[0,1,0]; 
color = [.7 .7 .7];
newnames = ["SP500", "EBP", "GZ"];

hor_back =1;
hor = hor-1;
fig = figure;
for k=1:n
    subplot(2,2,k)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k+n,:), fliplr(lowproj90(k+n,:))], color1, EdgeColor="none", FaceAlpha = 0.3)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k,:), fliplr(lowproj90(k,:))], color2, EdgeColor="none", FaceAlpha = 0.3)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k+n,:), fliplr(lowproj95(k+n,:))], color1, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(k+n,:), "MarkerFaceColor",color1)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k,:), fliplr(lowproj95(k,:))], color2, EdgeColor="none", FaceAlpha = 0.1)
hold on
plot(-hor_back:hor, proj(k,:), "MarkerFaceColor",color2)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(k),'interpreter', 'latex')
hold off
end
legend('$\Delta n_t <0$','$\Delta n_t \geq 0$','interpreter', 'latex')
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Months','interpreter', 'latex');
%print -depsc nonlin.eps

fig1 = figure;
    subplot(2,3,1)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(1,:), fliplr(lowproj90(1,:))], color, EdgeColor="none", FaceAlpha = 0.5)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(1,:), fliplr(lowproj95(1,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
plot(-hor_back:hor, proj(1,:), "MarkerFaceColor",color2)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(1),'interpreter', 'latex')
ylabel('$n_t \geq 0$','interpreter', 'latex');
yticks(linspace(-0.04,0.06,6))
hold off
for k=2:n
    subplot(2,3,k)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k,:), fliplr(lowproj90(k,:))], color, EdgeColor="none", FaceAlpha = 0.5)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k,:), fliplr(lowproj95(k,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
plot(-hor_back:hor, proj(k,:), "MarkerFaceColor",color2)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(k),'interpreter', 'latex')
hold off
end
    subplot(2,3,n+1)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(n+1,:), fliplr(lowproj90(n+1,:))], color, EdgeColor="none", FaceAlpha = 0.5)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(n+1,:), fliplr(lowproj95(n+1,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
plot(-hor_back:hor, proj(n+1,:), "MarkerFaceColor",color2)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(1),'interpreter', 'latex')
ylabel('$n_t < 0$','interpreter', 'latex');
hold off
for k=n+2:2*n
    subplot(2,3,k)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj90(k,:), fliplr(lowproj90(k,:))], color, EdgeColor="none", FaceAlpha = 0.5)
hold on
fill([-hor_back:hor, fliplr(-hor_back:hor)], [upproj95(k,:), fliplr(lowproj95(k,:))], color, EdgeColor="none", FaceAlpha = 0.3)
hold on
plot(-hor_back:hor, proj(k,:), "MarkerFaceColor",color2)
hold on
plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
box on
grid on
xlim([-hor_back,hor])
title(newnames(k-n),'interpreter', 'latex')
hold off
end
han=axes(fig1,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
%han.YLabel.Visible='on';
%ylabel(han,'p.p.','interpreter', 'latex');
xlabel(han,'Months','interpreter', 'latex');
%print -depsc nonlin_tiled_jlnm.eps