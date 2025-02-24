function U = expeuler_mat_syl(U,tau,m,A,B,g)

  E1 = expm(tau*A);
  E2 = expm(tau*B);
  for jj = 1:m
    gU = g(U);
    tmp = sylvester(tau*A,tau*B.',E1*(gU*E2.') - gU);
    U = E1*(U*E2.') + tau*tmp;
  end
end
