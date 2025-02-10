clear all
close all
T = 40;
lambda = [-1;-100];
A = diag(lambda);
y0 = [1;1];
FE_ts_restriction = min(-2*real(lambda)./abs(lambda).^2)
I = eye(2);
ts = 3000;
k = T/ts
yFE = y0;
yBE = y0;
yTRAP = y0;
for n = 1:ts
  yFE(:,n+1) = yFE(:,n)+k*A*yFE(:,n);
  yBE(:,n+1) = (I-k*A)\yBE(:,n);
  yTRAP(:,n+1) = (I-k/2*A)\((I+k/2*A)*yTRAP(:,n));
end
t = linspace(0,T,ts+1);
plot(t,max(abs(yFE)),'*',t,max(abs(yBE)),'o',t,max(abs(yTRAP)),'x')
xlabel('time')
ylabel('Solution infinity norm')
legend('Forward Euler','Backward Euler','Trapezoidal')
tsrange = 2100:100:ts;
counter = 0;
for ts = tsrange
  counter = counter+1;
  k = T/ts;
  yFE = y0;
  eFE(1) = 0;
  yBE = y0;
  eBE(1) = 0;
  yTRAP = y0;
  eTRAP(1) = 0;
  yEX = y0;
  for n = 1:ts
    yFE(:,n+1) = yFE(:,n)+k*A*yFE(:,n);
    yBE(:,n+1) = (I-k*A)\yBE(:,n);
    yTRAP(:,n+1) = (I-k/2*A)\((I+k/2*A)*yTRAP(:,n));
  end
  yEX = diag(exp(T*lambda))*y0;
  eFE(counter) = norm(yFE(:,ts+1)-yEX,inf)/norm(yEX,inf);
  eBE(counter) = norm(yBE(:,ts+1)-yEX,inf)/norm(yEX,inf);
  eTRAP(counter) = norm(yTRAP(:,ts+1)-yEX,inf)/norm(yEX,inf);
end
figure
loglog(tsrange,eFE,'*',tsrange,eBE,'o',tsrange,eTRAP,'x')
hold on
loglog([2000,2000],[1e-50,1e300])
xlabel('time')
ylabel('Error infinity norm')
legend('Forward Euler','Backward Euler','Trapezoidal')
