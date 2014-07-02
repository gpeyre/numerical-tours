function [h,t] = hist(x,n)

// hist - compute and display histogram
//
//  [h,t] = hist(x,n);
//
//  if h is ommitted, then display the histogram, 
//  otherwise, compute it.
//
//  Copyright (c) 2008 Gabriel Peyre


if argn(2)<2
    n = 100;
end

if argn(1)==0
    if length(n)>1
        n = length(n);
    end
    histplot(n,x(:));
    return;
end


x = x(:);
if length(n)>1
    t = n;
    n = length(n);
    // closest point histograms
    h = zeros(n,1);
    for i=1:length(x)
        [tmp,j] = compute_min( abs(x(i)-t(:)), 1 );
        h(j) = h(j)+1;
    end
else
    // equispaced histograms
    a = min(x); b = max(x);
    tau = (b-a)/n;
    a1 = a+tau/2; b1 = b-tau/2;
    t = a1:tau:b1;
    x1 = (x-a1)/(b1-a1)*(n-1)+1;
    x1 = round(x1);
    h = zeros(n,1);
    for i=1:n
        h(i) = sum(x1==i);
    end
end
h = h/sum(h);



endfunction
    

    
    


