function ret = phi2(z)
  idx = abs(z) < 1;
  ret = zeros (size (z));
  a = 1./factorial(17:-1:2);
  ret(idx) = polyval(a,z(idx));

  ret(~idx) = (exp(z(~idx)) - z(~idx) -1)./(z(~idx).^2);
