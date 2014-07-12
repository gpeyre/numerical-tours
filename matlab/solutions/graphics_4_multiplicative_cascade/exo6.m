sigma2list = [0.02 0.02 0.2 0.2];
mult = [1 10 1 1];
mulist = -mult.*sigma2list/2;
k = 0;
clf;
for k=1:4
   mu = mulist(k);
   sigma2 = sigma2list(k);
   if -(exp(2*(sigma2+mu))-1)<-1  % Condition of non-degeneracy
       disp('Be careful ! This cascade will degenerate as rmin -> 0 !')
   end
   H1 = 1 - exp(mu+sigma2/2);
   Wi = exp( randn(N,1)*sqrt(sigma2)+mu );
   f = ones(n)/rmin^H1;
   for i=1:N
       %progressbar(i,N);
       I = find( (X-xi(i)).^2+(Y-yi(i)).^2 <=ri(i)^2/4 );
       f(I) = f(I) * Wi(i);
   end
   subplot(2,2,k);
   imageplot(f, ['mu=' num2str(mu) ', sigma^2=' num2str(sigma2)]);
end
%EXO