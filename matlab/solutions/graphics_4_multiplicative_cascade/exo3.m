sigma2list = [0.1 0.1 0.5 0.5];
mult = [1 10 1 10];
mulist = -mult.*sigma2list;
ti = -1/2+rand(1,N) * (T+1);
ti_1 = -1/2+rand(1,round(T+1))*(T+1);
ti = [ti_1 ti];   % for exact scale invariance
umax = 1/rmin-1;
ui = [zeros(1,length(ti_1)) rand(1,N) * umax ]; % ui = 1-1/ri^{1+beta}
ri = (1+ui).^(-1);
clf
for k=1:4
   mu = mulist(k);
   sigma2 = sigma2list(k);
   if -(exp(2*(sigma2+mu))-1)<-1  % Condition of non-degeneracy
       disp('Be careful ! This cascade will degenerate as rmin -> 0 !')
   end
   H1 = 1 - exp(mu+sigma2/2);
   Wi = exp( randn(N,1)*sqrt(sigma2)+mu );
   f = ones(n,1)* exp(H1) / rmin^H1;
   for i=1:N
       I = find(abs(t-ti(i))<=ri(i)/2);
       f(I) = f(I) * Wi(i);
   end
   subplot(2,2,k);
   plot(t,f)
   axis([0 T 0 1.1*max(f)]) 
   title(['mu=' num2str(mu) ', sigma^2=' num2str(sigma2)]);
end
