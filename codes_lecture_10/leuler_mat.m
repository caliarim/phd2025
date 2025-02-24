function U = leuler_mat(U,tau,m,A,B,g)

  E1 = expm(tau*A);
  E2 = expm(tau*B);
  for jj = 1:m
    U = E1*((U+tau*g(U))*E2.');
  end
end
