cd 'C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty'
addpath("C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\Data")
addpath("C:\Users\giaco\Desktop\phd\Useful data")

close all
clear
data = readtable("Levels.xlsx");
X = table2array(data(:,2:end));
tobefiltered = X(:,[1,2,3,4,6,19]);
T = zeros(size(tobefiltered));
C = zeros(size(tobefiltered));
for i=1:size(tobefiltered,2)
    [T(:,i), C(:,i)] = hpfilter(tobefiltered(:,i));
end
[Tinv,Cinv] = hpfilter(X(:,5)); %GFCF
X(:,[1,2,3,4,6,19]) = C(:,:)./T(:,:)*100;
X(5:end,5) = Cinv./Tinv*100;

for i =1:4
X(:,i) = movmean(X(:,i),2);
end

controls = X(:,[17:19]);
hor =12;
p=0;
hor_back = 0;
n = 4;
proj = zeros(2*n,hor+hor_back+1);
upproj95 = zeros(2*n,hor+hor_back+1);
lowproj95 = zeros(2*n,hor+hor_back+1);
upproj90 = zeros(2*n,hor+hor_back+1);
lowproj90 = zeros(2*n,hor+hor_back+1);

shock = [std(X(:,end-3)), std(X(:,end-1)), std(X(:,end))];
capratio = lagmatrix(X(:,9),1);% = double(X(:,9)<0);
z       = rmmissing([X(:,[1:4,10,11,12,16,19]),X(:,9), X(:,end-3), X(:,9).*X(:,end-3)]);
zlags   = lagmatrix(z(:, 1:end-2), 0:p);
xx      = [zlags(p+1:end,:), z(p+1:end, end-1:end)];
r_mean = mean(capratio(capratio>=0));%mean(rmmissing(controls(:,end)));
rr = mean(capratio(capratio<0));%mean(rmmissing(controls(:,end))) - std(rmmissing(controls(:,end)));


for j = 1:n
    for h=1:hor
        %lags_shock = rmmissing(lagmatrix(xx(:,end-2:end-1),1:h));
        %zh =  [lags_shock(1:end,end-3), xx(1:end-h,2:end)];
        zh = [xx(1:end-h,1:end)];
        %zh(:,j) = [];
        %zh = [xx(h:end-1,2:end-3), xx(1:end-h,end-2:end),xx(1:end-h,end).*xx(1:end-h,11)];
        x = xx(h+1:end,j);
        lp = fitlm(zh,x);
       aa = shock(1)*table2array(lp.Coefficients(end,1));
       
        bb      =               shock(1)*table2array(lp.Coefficients(end-1,1));
       
        proj(j,h+hor_back+1)      =     bb+aa*r_mean;
        varproj = r_mean^2*table2array(lp.Coefficients(end,2))^2 + table2array(lp.Coefficients(end-1,2))^2 ...
                +2*lp.CoefficientCovariance(end,end-1)*r_mean;
        upproj95(j,h+hor_back+1) = proj(j,h+1) ...
            + shock(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j,h+hor_back+1) = proj(j,h+1) ...
            - shock(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j,h+hor_back+1) = proj(j,h+1)...
            + shock(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j,h+hor_back+1) = proj(j,h+1)...
            - shock(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));

        proj(j+n,h+hor_back+1)      =     bb+aa*rr;
        varproj = rr^2*table2array(lp.Coefficients(end,2))^2 + table2array(lp.Coefficients(end-1,2))^2 ...
                +2*lp.CoefficientCovariance(end,end-1)*rr;
        upproj95(j+n,h+hor_back+1) = proj(j+n,h+1) ...
            + shock(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j+n,h+hor_back+1) = proj(j+n,h+1) ...
            - shock(1)*varproj^(1/2)*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j+n,h+hor_back+1) = proj(j+n,h+1)...
            + shock(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j+n,h+hor_back+1) = proj(j+n,h+1)...
            - shock(1)*varproj^(1/2)*tinv(0.95,size(zh,1)-size(zh,2));
    end
end
color1 =  [0.4940 0.1840 0.5560];%[0 0.4470 0.7410];%[0,0,1];
color2 = [0.3010 0.7450 0.9330];%[0,1,0]; 
color = [.7 .7 .7];

newnames= ["Credit", "Treasuries","Net Worth", "Debt Instruments"];
%newnames= ["Investment", "Consumption", "Hours Worked"];
hor_back=1;
hor=hor-1;
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
xlabel(han,'Quarters','interpreter', 'latex');
%print -depsc nonlin_realsector.eps

fig1 = figure;
    subplot(2,n,1)
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
hold off
for k=2:n
    subplot(2,n,k)
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
    subplot(2,n,n+1)
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
    subplot(2,n,k)
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
xlabel(han,'Quarters','interpreter', 'latex');
%print -depsc nonlin_finsector_tiled_jlnm.eps