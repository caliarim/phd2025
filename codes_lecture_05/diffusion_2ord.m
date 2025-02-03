clear all
close all

N = 200;

x = linspace(0,1,N+2)';
x = x(2:N+1); % inner nodes

h = 1/(N+1);

y0 = 4*x.*(1-x);

delta = 1;

D2 = toeplitz([-2,1,zeros(1,N-2)]/(h^2));

A = delta*D2;

[V,D] = eig(A);

d = diag(D);

g = @(y) 1./(1+y.^2); % nonlinear part

tstar = 0.1;

mref = 10000;
m = mref;
tau = tstar/m;

yref = erk_p1(1/2,y0,m,tau,g,V,d);
normref = norm(yref,inf);

figure;
plot(x,y0,'xb',x,yref,'or')
legend('Initial','Final')

mrange = 10:10:100;
counter = 0;

c2_p1 = 1/2;

c2_p2 = 1;

for m = mrange
  counter = counter + 1;
  tau = tstar/m;

  y_erk_p1 = erk_p1(c2_p1,y0,m,tau,g,V,d);
  err_erk_p1(counter) = norm(yref-y_erk_p1,inf)/normref;

  y_erk_p1p2 = erk_p1p2(c2_p2,y0,m,tau,g,V,d);
  err_erk_p1p2(counter) = norm(yref-y_erk_p1p2,inf)/normref;

end

figure;
loglog(mrange,err_erk_p1,'xb')
hold on
loglog(mrange,err_erk_p1p2,'or')
loglog(mrange, err_erk_p1(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange, err_erk_p1(end)*(mrange/mrange(end)).^(-2),'-k')
loglog(mrange, err_erk_p1p2(end)*(mrange/mrange(end)).^(-2),'-k')
legend('Expint phi1','Expint phi1phi2')
