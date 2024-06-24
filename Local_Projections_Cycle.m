
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

controls = X(:,[17:19]);
hor =12;
p=2;
hor_back = 0;
n = 4;
proj = zeros(3*n,hor+hor_back+1);
upproj95 = zeros(3*n,hor+hor_back+1);
lowproj95 = zeros(3*n,hor+hor_back+1);
upproj90 = zeros(3*n,hor+hor_back+1);
lowproj90 = zeros(3*n,hor+hor_back+1);

shock = [std(X(:,end-3)), std(X(:,end-1)), std(X(:,end))];
%X(:,[end-3,end-1:end]) =X(:,[end-3,end-1:end])/std(X(:))
z       = rmmissing([X(:,[1:4,10,11,12,16,19]), X(:,[end-3,end-1:end])]);
zlags   = lagmatrix(z(:, 1:end-3), 0:p);
xx      = [zlags(p+1:end,:), z(p+1:end, end-2:end)];
for i=0:2
for j = 1:n
    for h=1:hor
        lags_shock = rmmissing(lagmatrix(xx(:,end-2:end), 1:h));
        zh = xx(1:end-h,[1:end-3,end-i]);
        %zh = [lags_shock(:,1:end-3), xx(1:end-h,:)];
        %zh = [xx(h:end-1,1:end-3), xx(1:end-h,end)];
        %zh(:,j) = [];
        x = z(p+1+h:end,j);
        lp = fitlm(zh,x);
        proj(j+i*n,h+hor_back+1) = shock(end-i)*table2array(lp.Coefficients(end,1));
        upproj95(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) + shock(end-i)*table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        lowproj95(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) -shock(end-i)* table2array(lp.Coefficients(end,2))*tinv(0.975,size(zh,1)-size(zh,2));
        upproj90(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) + shock(end-i)*table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));
        lowproj90(j+i*n,h+hor_back+1) = proj(j+i*n,h+1) -shock(end-i)* table2array(lp.Coefficients(end,2))*tinv(0.95,size(zh,1)-size(zh,2));

    end
end
end

newnames= ["Credit", "Treasuries", "Net Worth"];
%newnames= ["Investment", "Consumption", "Hours Worked"];
%fig = figure;
%for k=1:n
%    subplot(1,3,k)
%hold on
%fill([-hor_back:hor, fliplr(-hor_back:hor)], [upprojf90(k,:), fliplr(lowprojf90(k,:))], [0.3010 0.7450 0.9330], EdgeColor="none", FaceAlpha = 0.3)
%hold on
%fill([-hor_back:hor, fliplr(-hor_back:hor)], [upprojm90(k,:), fliplr(lowprojm90(k,:))], [0.8500 0.3250 0.0980], EdgeColor="none", FaceAlpha = 0.3)
%hold on
%fill([-hor_back:hor, fliplr(-hor_back:hor)], [upprojf95(k,:), fliplr(lowprojf95(k,:))], [0.3010 0.7450 0.9330], EdgeColor="none", FaceAlpha = 0.1)
%hold on
%plot(-hor_back:hor, proj_f(k,:), "MarkerFaceColor",[0.8500 0.3250 0.0980])
%hold on
%fill([-hor_back:hor, fliplr(-hor_back:hor)], [upprojm95(k,:), fliplr(lowprojm95(k,:))], [0.8500 0.3250 0.0980], EdgeColor="none", FaceAlpha = 0.1)
%hold on
%plot(-hor_back:hor, proj_m(k,:), "MarkerFaceColor",[0.3010 0.7450 0.9330])
%hold on
%plot(-hor_back:hor, zeros(1,hor_back+hor+1), 'r')
%box on
%grid on
%xlim([-hor_back,hor])
%title(newnames(k),'interpreter', 'latex')
%hold off
%end
%legend('Financial','Macro')
%han=axes(fig,'visible','off'); 
%han.Title.Visible='on';
%han.XLabel.Visible='on';
%han.YLabel.Visible='on';
%ylabel(han,'p.p.','interpreter', 'latex');
%xlabel(han,'Quarters','interpreter', 'latex');
%print -depsc lv_realsector.eps
color = [.7 .7 .7];
hor_back = 1;
hor=hor-1;
fig1 = figure;
tcl = tiledlayout(3,3);
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
han.Title.Visible='on';
%han.XLabel.Visible='on';
%han.YLabel.Visible='on';
%ylabel(han,'p.p.','interpreter', 'latex');
%xlabel(han,'Quarters','interpreter', 'latex');
%print -depsc lv_finsector_tiled.eps

