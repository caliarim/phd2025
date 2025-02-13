function y = exp_taylor(z)
  tol = eps;
  y = 1;
  w = 1;
  m = 0;
  while (abs(w) > tol * y)
    m = m+1;
    w = w*z/m;
    y = y+w
  end
