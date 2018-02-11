f = ones(n,1) * exp(H1) / rmin^H1;
k = 0;
clf
for i=1:N
   I = find(abs(t-ti(i))<=ri(i)/2);
   f(I) = f(I) * Wi(i);    
   if i==10 || i==80 || i==round(.1*N) || i==round(.5*N)
       k = k+1;
       subplot(2,2,k);
       plot(t,f)
		axis([0 T 0 1.1*max(f)])
   end
end
