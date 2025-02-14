function [a,b] = phi_coeff(ell,p)
%
% function [a,b] = phi_coeff(ell,p)
%
% Returns the Pade' coefficients
%
% a(1)+a(2)*z+...+a(p)*z^(p-1)+a(p+1)*z^p
% ---------------------------------------
% b(1)+b(2)*z+...+b(p)*z^(p-1)+b(p+1)*z^p
%
% of the rational approximation of phi_ell
ell = sym(ell);
for i = 0:p
  j = 0:i;
  a(i+1)= factorial(p)/factorial(2*p+ell)*sum(factorial(2*p+ell-j).*(-1).^j./(factorial(j).*factorial(p-j).*factorial(ell+i-j)));
  b(i+1) = factorial(p)/factorial(2*p+ell)*factorial(2*p+ell-i)/(factorial(i)*factorial(p-i))*(-1)^i;
end
