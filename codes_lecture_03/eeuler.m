function y = eeuler(m,tau,A,f,y0)

y = y0;
for jj=1:m
  y = y + tau*f(y);
end
