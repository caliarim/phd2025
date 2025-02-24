clear all
close all

% 2D ADR with homogeneous Dirichlet/Neumann b.c.

nx = 50;
ny = 51;

delta = 0.05;
alpha = -0.8;

mref = 5000;

mrange = 100:10:200;

T = 0.6;

g = @(u) 1./(1+u.^2);

x = linspace(0,1,nx+1).';
x = x(2:nx+1);
y = linspace(0,1,ny+1).';
y = y(2:ny+1);

hx = 1/nx;
hy = 1/ny;

Pe = max(abs(alpha)*[hx,hy]/(2*delta))

[X,Y] = ndgrid(x,y);

Ix = speye(nx);
Iy = speye(ny);

D1x = spdiags(ones(nx,1)*[-1,0,1]/(2*hx),-1:1,nx,nx);
D1x(nx,nx-1) = 0;
D2x = spdiags(ones(nx,1)*[1,-2,1]/(hx^2),-1:1,nx,nx);
D2x(nx,nx-1) = 2/(hx^2);
D1y = spdiags(ones(ny,1)*[-1,0,1]/(2*hy),-1:1,ny,ny);
D1y(ny,ny-1) = 0;
D2y = spdiags(ones(ny,1)*[1,-2,1]/(hy^2),-1:1,ny,ny);
D2y(ny,ny-1) = 2/(hy^2);

A = delta*D2x + alpha*D1x;
B = delta*D2y + alpha*D1y;

K = kron(Iy,A) + kron(B,Ix);

%[V,D] = eig(full(A));
%cond(V)
%pause

U0 = 256*(X.*(1-X).*Y.*(1-Y)).^2;

figure;
mesh(X,Y,U0)
title('Initial')
xlabel('x'); ylabel('y'); zlabel('u(x,y)')
drawnow

uref = leuler_mat(U0,T/mref,mref,A,B,g);
normref = norm(uref(:),inf);

figure;
mesh(X,Y,reshape(uref,nx,ny));
title('Final')
xlabel('x'); ylabel('y'); zlabel('u(x,y)')
drawnow

counter = 0;
for m = mrange
  counter = counter + 1;
  m
  tau = T/m;
  u_lev = leuler_vec(U0(:),tau,m,K,g);
  err_lev(counter) = norm(uref(:)-u_lev,inf);
  u_lem = leuler_mat(U0,tau,m,full(A),full(B),g);
  err_lem(counter) = norm(uref(:)-u_lem(:),inf);
  u_eev = expeuler_vec(U0(:),tau,m,K,g);
  err_eev(counter) = norm(uref(:)-u_eev,inf);
  u_eems = expeuler_mat_syl(U0,tau,m,full(A),full(B),g);
  err_eems(counter) = norm(uref(:)-u_eems(:),inf);
  u_eemsp = expeuler_mat_split(U0,tau,m,full(A),full(B),g);
  err_eemsp(counter) = norm(uref(:)-u_eemsp(:),inf);

end

figure;
loglog(mrange,err_lev,'xb',mrange,err_lev(end)*(mrange/mrange(end)).^(-1),'--k')
hold on
loglog(mrange,err_lem,'sm',mrange,err_lem(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange,err_eev,'or',mrange,err_eev(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange,err_eems,'+g',mrange,err_eems(end)*(mrange/mrange(end)).^(-1),'--k')
loglog(mrange,err_eemsp,'*c',mrange,err_eemsp(end)*(mrange/mrange(end)).^(-1),'--k')
drawnow
