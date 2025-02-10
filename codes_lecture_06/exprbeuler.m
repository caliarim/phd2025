function y = exprbeuler(y0,m,tau,f,JJ)

  y = y0;
  for jj = 1:m
    [V,D] = eig(JJ(y));
    P1 = diag(phi1(tau*diag(D)));
    y = y + tau*(V*(P1*(V'*f(y))));
  end

end
