clear all
close all
m = [30,31];
N = prod(m);
x = linspace(0,1,m(1)+2)';
x = x(2:m(1)+1);
h(1) = 1/(m(1)+1);
y = linspace(0,1,m(2)+2)';
y = y(2:m(2)+1);
h(2) = 1/(m(2)+1);
A1 = spdiags(ones(m(1),1)*[1,-2,1]/h(1)^2,-1:1,m(1),m(1));
A2 = spdiags(ones(m(2),1)*[1,-2,1]/h(2)^2,-1:1,m(2),m(2));
A = kron(speye(m(2)),A1)+kron(A2,speye(m(1)));
[X,Y] = ndgrid(x,y);
U0 = sin(pi*X).*sin(2*pi*Y);
u0 = U0(:);
mesh(X,Y,reshape(A*u0,size(X)))
xlabel('x')
ylabel('y')
title('A*u0')
figure
mesh(X,Y,-U0*pi^2-U0*4*pi^2)
xlabel('x')
ylabel('y')
title('Laplacian of u0')
tau = 0.1;
tic
r1 = expm(tau*A)*u0;
toc
tic
r2 = expmv(A,u0,tau); % https://github.com/higham/expmv
toc
tic
r3 = kiops(tau,A,u0,[],[],[],[],false); % https://gitlab.com/stephane.gaudreault/kiops
toc
norm(r1-r2,inf)
norm(r1-r3,inf)
u1 = u0/2;
u2 = u0/4;
Ahat = [A,u2/tau,u1;spalloc(2,N,0),spdiags([1;1],1,2,2)];
tic
r1 = expm(tau*A)*u0+tau*phi1m(tau*A)*u1+tau*phi2m(tau*A)*u2;
toc
tic
r2 = expmv(Ahat,[u0;0;1],tau);
toc
r2 = r2(1:N);
tic
r3 = kiops(tau,A,[u0,u1,u2/tau],[],[],[],[],false);
toc
norm(r1-r2,inf)
norm(r1-r3,inf)

