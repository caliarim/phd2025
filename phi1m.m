function [N,PHI0] = phi1m(A)
%
% function PHI1 = phi1m(A)
% function [PHI1,PHI0] = phi1m(A)
%
% N is phi_1(A) and PHI0 is equivalent to expm(A)

% Scale A by power of 2 so that its norm is < 1/2 .
[f,e] = log2(norm(A,1));
s = min(max(0,e+1),1023);
A = A/2^s;
% Pade approximation for phi1(z)
ID = eye(size(A));
a = [1/259459200,1/3603600,1/171600,1/4680,1/936,1/30,1/30,1];
b = [-1/32432400,1/514800,-1/17160,1/936,-1/78,1/10,-7/15,1];
p = length(a);
% Horner's scheme
N = a(1)*A+a(2)*ID;
D = b(1)*A+b(2)*ID;
for i = 3:p
   N = N*A+a(i)*ID;
   D = D*A+b(i)*ID;
end
N = full(D\N);
% Undo scaling by repeated squaring
PHI0 = A*N+ID;
for i = 1:s
  N = (PHI0+ID)*N/2;
  PHI0 = PHI0*PHI0;
end
