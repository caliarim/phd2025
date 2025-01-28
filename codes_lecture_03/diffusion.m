clear all
close all

N = 100;

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
f = @(y) A*y + g(y); % rhs

tstar = 0.1;

mref = 10000;
m = mref;
tau = tstar/m;

yref = expeuler(y0,m,tau,f,V,d);
normref = norm(yref,inf);

figure;
plot(x,y0,'xb',x,yref,'or')
legend('Initial','Final')

mrange = 1000:100:5000;
counter = 0;

for m = mrange
  counter = counter + 1;
  tau = tstar/m;

  y_expeuler = expeuler(y0,m,tau,f,V,d);
  err_expeuler(counter) = norm(yref-y_expeuler,inf)/normref;
  y_bfeuler = bfeuler(y0,m,tau,g,V,d);
  err_bfeuler(counter) = norm(yref-y_bfeuler,inf)/normref;

  y_ee = ee(y0,m,tau,f);
  err_ee(counter) = norm(yref-y_ee,inf)/normref;

end

figure;
loglog(mrange,err_expeuler,'xb')
hold on
loglog(mrange,err_bfeuler,'or')
loglog(mrange,err_ee,'sg')
loglog(mrange, err_expeuler(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange, err_bfeuler(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange, err_ee(end)*(mrange/mrange(end)).^(-1),'--k')
legend('Exp Euler','BF Euler','Explicit Euler')
