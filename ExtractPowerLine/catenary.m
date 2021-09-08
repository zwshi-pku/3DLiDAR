%% Catenary: y=a+ccosh((x-b)/c)
function cfun = catenary(xs,ys)
syms x;
f = fittype('a+c*cosh((x-b)/c)','independent','x','coefficients',{'a','c','b'});
cfun = fit(xs,ys,f,'StartPoint', [-780.6,780.4,-1.184])
end