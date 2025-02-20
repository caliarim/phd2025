function [V,H,vm1] = arnoldi(A,v,m)
%
% Arnoldi factorization by Modified Gram-Schmidt
%
% V'*A*V = H(1:m,1:m)
% A*V = [V,vm1]*H
%
% m must be smaller or equal than length(A)
  if (m < 1)
    error('arnoldi: m must be larger than 0')
  end
  H = zeros(m+1,m);
  V = zeros(size(A,1),m+1);
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
  else
    vm1 = V(:,m+1);
  end
  V = V(:,1:m);
%!test
%! A = randn(10);
%! v = randn(10,1);
%! m = 5;
%! [V,H,vm1] = arnoldi(A,v,m);
%! assert(V'*A*V,H(1:m,1:m),8*eps)
%! assert(A*V,[V,vm1]*H,8*eps)
%! assert(V'*V,eye(m),2*eps)
%!test
%! A = randn(10);
%! v = randn(10,1);
%! m = 10;
%! [V,H,vm1] = arnoldi(A,v,m);
%! assert(V'*A*V,H(1:m,1:m),256*eps)
%! assert(A*V,[V,vm1]*H,16*eps)
%! assert(V'*V,eye(m),16*eps)
