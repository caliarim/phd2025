clear all
close all

N = 100;

xl = 0;
xr = 1;

delta = 1;

T = 0.1;
mrange = 10:10:100;
mref = 10000;

% Homogeneous Dirichlet
x = linspace(xl,xr,N+2).';
x = x(2:N+1);
h = (xr-xl)/(N+1);

y0 = 4*x.*(1-x);

figure
plot(x,y0,'xb')
title('Initial condition')

D2 = toeplitz([-2,1,zeros(1,N-2)]/(h^2));

A = delta*D2;

[V,D] = eig(A); % A*V = V*D

disp(sprintf('The condition number is %.3f',cond(V)))

%r1 = expm(tau*A)*y0;
%r2 = V*(diag(exp(diag(tau*D)))*(V\y0));
%r2_b = V*(diag(exp(diag(tau*D)))*(V'*y0));
%norm(r1-r2,inf)/norm(r1,inf)

g = @(y) 1./(1+y.^2);
f = @(y) A*y + g(y);

%yref = expeuler(mref,T/mref,A,f,y0,V,D);
yref = bfeuler(mref,T/mref,A,g,y0,V,D);
normref = norm(yref,inf);
figure
plot(x,y0,'xb')
hold on
plot(x,yref,'or')
legend('Initial solution','Final solution')

counter = 0;

for m = mrange
  counter = counter + 1;
  tau = T/m;
  y_bfe = bfeuler(m,tau,A,g,y0,V,D);
  y_expe = expeuler(m,tau,A,f,y0,V,D);
  %y_ee = eeuler(m,tau,A,f,y0);
  err_bfe(counter) = norm(y_bfe-yref,inf)/normref;
  err_expe(counter) = norm(y_expe-yref,inf)/normref;
  %err_ee(counter) = norm(y_ee-yref,inf)/normref;
end

figure;
loglog(mrange,err_bfe,'xm')
hold on
loglog(mrange,err_expe,'+g')
%loglog(mrange,err_ee,'^b')
loglog(mrange,err_bfe(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange,err_expe(end)*(mrange/mrange(end)).^(-1),'--k')
%loglog(mrange,err_ee(end)*(mrange/mrange(end)).^(-1),'--k')
legend('BF euler','EXP euler','E euler')
