function u = leuler_vec(u,tau,m,K,g)

  tol = tau^2/100;
  mkry = 10;
  for jj = 1:m
    [u,mkry] = kiops(tau,K,u+tau*g(u),tol,mkry,10,128,false);
  end
end
