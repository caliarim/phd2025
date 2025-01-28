function y = bfeuler(y0,m,tau,g,V,d)

  y = y0;
  M = diag(1./(1-tau*d));
  for jj = 1:m
    y = V*(M*(V'*(y + tau*g(y))));
  end

end
