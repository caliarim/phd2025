clear all
close all

N = 100;
tstar = 0.1;
delta = 1i;

x = linspace(0,1,N+2)';
x = x(2:N+1); % inner nodes
h = 1/(N+1);

y_ex = @(t) exp(t)*(4*x.*(1-x));
%y_ex = @(t) exp(t)*sin(2*pi*x);

g = @(t) 4*exp(t)*(x.*(1-x)+2*delta);
%g = @(t) exp(t)*(1+delta*(2*pi)^2)*sin(2*pi*x);

y0 = y_ex(0);

D2 = toeplitz([-2,1,zeros(1,N-2)]/(h^2));

A = delta*D2;

[V,D] = eig(A);

d = diag(D);

ysol = y_ex(tstar);
%mref = 10000;
%ysol = expquad2(0,1,y0,mref,tstar/mref,g,A,V,d);
normsol = norm(ysol,inf);

figure;
plot(x,y0,'xb',x,ysol.*conj(ysol),'or')
legend('Initial','Final')

mrange = 1:1:200;
counter = 0;

% SINGLE POINT
c1s = 0;

% TRAPEZOIDAL
%c1 = 0;
%c2 = 1;

% GAUSS RADAU
%c1 = 0;
%c2 = 2/3;

% GAUSS
%c1 = (-1/sqrt(3)+1)/2;
%c2 = (1/sqrt(3)+1)/2;

for m = mrange
  m
  counter = counter + 1;
  tau = tstar/m;

  y_eq1 = expquad1(c1s,y0,m,tau,g,A,V,d);
  err_eq1(counter) = norm(ysol-y_eq1,inf)/normsol;

  %y_eq2 = expquad2(c1,c2,y0,m,tau,g,A,V,d);
  %err_eq2(counter) = norm(ysol-y_eq2,inf)/normsol;

end

figure;
loglog(mrange,err_eq1,'xb')
hold on
%loglog(mrange,err_eq2,'or')
loglog(mrange, err_eq1(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange, err_eq1(end)*(mrange/mrange(end)).^(-2),'-k')
%loglog(mrange, err_eq2(1)*(mrange/mrange(1)).^(-3),'-.k')
%legend('ExpQuad1','ExpQuad2')
