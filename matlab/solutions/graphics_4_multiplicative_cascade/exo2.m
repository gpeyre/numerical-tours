clf
for k=1:4
   Wi = exp( randn(N,1)*sqrt(sigma2)+mu );
   ti = -1/2+rand(1,N) * (T+1);
   ti_1 = -1/2+rand(1,round(T+1))*(T+1);
   ti = [ti_1 ti];   % for exact scale invariance
   umax = 1/rmin-1;
   ui = [zeros(1,length(ti_1)) rand(1,N) * umax ]; % ui = 1/ri-1
   ri = (1+ui).^(-1);
   f = ones(n,1)* exp(H1) / rmin^H1;
   for i=1:N
       I = find(abs(t-ti(i))<=ri(i)/2);
       f(I) = f(I) * Wi(i);
   end
   subplot(2,2,k);
   plot(t,f)
   axis([0 T 0 1.1*max(f)])
end
