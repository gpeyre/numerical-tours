tv  = sum(sum( sqrt(sum( grad(f0).^2,3 )) ));
% display
clf;
imageplot(f0, strcat(['TV=', num2str(tv,4)]), 1,2,1);
