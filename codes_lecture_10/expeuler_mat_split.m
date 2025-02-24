function U = expeuler_mat_split(U,tau,m,A,B,g)

  [P1,E1] = phi1m(tau*A);
  [P2,E2] = phi1m(tau*B);
  for jj = 1:m
    gU = g(U);
    U = E1*(U*E2.') + tau*P1*(g(U)*P2.');
  end
end
