function y = erk_p1p2(c2,y0,m,tau,g,V,d)

  y = y0;
  Ec2 = diag(exp(tau*c2*d));
  P1c2 = diag(phi1(tau*c2*d));
  E = diag(exp(tau*d));
  P1 = diag(phi1(tau*d));
  P2 = diag(phi2(tau*d));
  for jj = 1:m
    y2 = V*(Ec2*(V'*y)) + c2*tau*(V*(P1c2*(V'*g(y))));
    y = V*(E*(V'*y)) + tau*(V*((P1-P2/c2)*(V'*g(y)))+V*((P2/c2)*(V'*g(y2))));
  end

end
