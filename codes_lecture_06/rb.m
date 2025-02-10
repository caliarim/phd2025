function y = rb(y0,m,tau,f,JJ)

  y = y0;
  for jj = 1:m
    [V,D] = eig(JJ(y));
    I1 = diag(1./(1-tau/2*diag(D)));
    k1 = tau*(V*(I1*(V'*f(y))));
    y = y + k1;
  end

end
