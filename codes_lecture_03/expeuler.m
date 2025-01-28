function y = expeuler(y0,m,tau,f,V,d)

  y = y0;
  M = diag(phi1(tau*d));
  for jj = 1:m
    y = y + tau*(V*(M*(V'*f(y))));
  end

end
