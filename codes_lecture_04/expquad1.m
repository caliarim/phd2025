function y = expquad1(c1,y0,m,tau,g,A,V,d)

  y = y0;
  M = diag(phi1(tau*d));
  t = 0;
  for jj = 1:m
    y = y + tau*(V*(M*(V'*(A*y+g(t+c1*tau)))));
    t = t + tau;
  end

end
