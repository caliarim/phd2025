function y = expquad2(c1,c2,y0,m,tau,g,A,V,d)

  y = y0;

  v1 = phi1(tau*d);
  v2 = phi2(tau*d);

  M0 = diag(exp(tau*d));
  M1 = diag(c2/(c2-c1)*v1-1/(c2-c1)*v2);
  M2 = diag(-c1/(c2-c1)*v1+1/(c2-c1)*v2);
  t = 0;
  for jj = 1:m
    y = V*(M0*(V'*y)) + tau*(V*(M1*(V'*g(t+c1*tau)))+V*(M2*(V'*g(t+c2*tau))));
    t = t + tau;
  end

end
