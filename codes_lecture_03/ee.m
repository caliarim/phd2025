function y = ee(y0,m,tau,f)
  y = y0;
  for jj = 1:m
    y = y + tau*f(y);
  end
end
