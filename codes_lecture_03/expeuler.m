function y = expeuler(m,tau,A,f,y0,V,D)

y = y0;
M = diag(phi1(diag(tau*D)));
for jj=1:m
  y = y + tau*(V*(M*V'*f(y)));
end
