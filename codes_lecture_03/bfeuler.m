function y = bfeuler(m,tau,A,g,y0,V,D)

y = y0;
M = diag(1./(1-tau*diag(D)));
for jj=1:m
  y = V*(M*V'*(y + tau*g(y)));
end
