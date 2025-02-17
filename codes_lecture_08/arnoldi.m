function [V,H,vm1] = arnoldi(A,v,m)
%
% Arnoldi factorization by Modified Gram-Schmidt
%
% V'*A*V = H(1:m,1:m)
% A*V = [V,vm1]*H
%
% m must be smaller or equal than length(A)
  V(:,1) = v/norm(v);
  for j = 1:m
    w = A*V(:,j);
    for i = 1:j
      H(i,j) = w'*V(:,i);
      w = w-H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    V(:,j+1) = w/H(j+1,j); % otherwise "breakdown"
  end
  if (m == length(A))
    warning('arnoldi: full factorization')
    vm1 = zeros(size(v));
    V = V(:,1:m);
  else
    vm1 = V(:,m+1);
    V = V(:,1:m);
  end
%!test
%! A = randn(10);
%! v = randn(10,1);
%! m = 5;
%! [V,H,vm1] = arnoldi(A,v,m);
%! assert(V'*A*V,H(1:m,1:m),8*eps)
%! assert(A*V,[V,vm1]*H,8*eps)
