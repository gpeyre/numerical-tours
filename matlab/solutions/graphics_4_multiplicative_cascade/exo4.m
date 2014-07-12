f = ones(n)/rmin^H1;
k = 0;
clf
for i=1:N  
   %progressbar(i,N);
   I = find( (X-xi(i)).^2+(Y-yi(i)).^2 <=ri(i)^2/4 );
   f(I) = f(I) * Wi(i); 
   if i==30 || i==200 || i==round(N*.5) || i==round(N*.8)  
       k = k+1;
       subplot(2,2,k);
       imageplot(f);
   end
end
