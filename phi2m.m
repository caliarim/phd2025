function [N,PHI1,PHI0] = phi2m(A)
%
% function PHI2 = phi2m(A)
% function [PHI2,PHI1] = phi2m(A)
% function [PHI2,PHI1,PHI0] = phi2m(A)

% Scale A by power of 2 so that its norm is < 1/2 .
[f,e] = log2(norm(A,1));
s = min(max(0,e+1),1023);
A = A/2^s;
% Pade approximation for phik(A)
ID = eye(size(A));
a = [1/121080960,0,3/80080,-1/2730,25/2184,-1/21,1/2];
b = [1/2162160,-1/40040,5/8008,-5/546,15/182,-3/7,1];
p = length(a);
N = a(1)*A+a(2)*ID;
D = b(1)*A+b(2)*ID;
for i = 3:p
   N = N*A+a(i)*ID;
   D = D*A+b(i)*ID;
end
N = full(D\N);
% Undo scaling by repeated squaring
PHI1 = A*N+ID;
PHI0 = A*PHI1+ID;
for i = 1:s
  N = ((PHI0+ID)*N+PHI1)/4;
  PHI1 = (PHI0+ID)*PHI1/2;
  PHI0 = PHI0*PHI0;
end
