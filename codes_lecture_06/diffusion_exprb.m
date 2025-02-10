clear all
close all

N = 150;

x = linspace(0,1,N+2)';
x = x(2:N+1); % inner nodes

h = 1/(N+1);

y0 = 4*x.*(1-x);

delta = 1;

D2 = toeplitz([-2,1,zeros(1,N-2)]/(h^2));

A = delta*D2;

g = @(y) 1./(1+y.^2); % nonlinear part

f = @(y) A*y + g(y);

JJ = @(y) A + diag(-2*y./((1+y.^2).^2));

tstar = 0.1;

mref = 2000;
m = mref;
tau = tstar/m;

yref = exprbeuler(y0,m,tau,f,JJ);
normref = norm(yref,inf);

figure;
plot(x,y0,'xb',x,yref,'or')
legend('Initial','Final')

mrange = 100:50:400;
counter = 0;

for m = mrange
  counter = counter + 1;
  tau = tstar/m;

  y_erb = exprbeuler(y0,m,tau,f,JJ);
  err_erb(counter) = norm(yref-y_erb,inf)/normref;

  y_rb = rb(y0,m,tau,f,JJ);
  err_rb(counter) = norm(yref-y_rb,inf)/normref;

end

figure;
loglog(mrange,err_erb,'xb')
hold on
loglog(mrange,err_rb,'or')
loglog(mrange, err_erb(end)*(mrange/mrange(end)).^(-2),'-k')
loglog(mrange, err_rb(end)*(mrange/mrange(end)).^(-2),'-k')
legend('Exp RB Euler', 'RB')
