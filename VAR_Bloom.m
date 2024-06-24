data = readtable("Levels.xlsx");

cd 'C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty'
addpath("C:\Users\giaco\Desktop\phd\Current Projects\Uncertainty\Data")
addpath("C:\Users\giaco\Desktop\phd\Useful data")

close all
clear
data = readtable("VARDATA.xlsx");
jlnshocks = readtable("unc_shock_monthly.csv");
X = table2array(data(:,3:end));
x=X(:,[7,8,9]);
x = rmmissing(x);
p =2;
xx = rmmissing(lagmatrix(x,1:p));
res= zeros(size(x,1)-p, size(x,2));
n = size(x, 2);
for i =1:n
    var = fitlm(xx, x(1+p:end,i));
    res(:,i) = x(1+p:end,i)-var.Fitted;
end

GG = res'*res/length(res);
G = chol(GG(1:n,1:n));
shocks = (G\res')';

jln = table2array(jlnshocks(19+p:570,2:end));
corrcoef(jln(:,end), shocks(:,2))