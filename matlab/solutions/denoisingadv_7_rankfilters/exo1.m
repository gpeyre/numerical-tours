beta_list = linspace(0,1, 6);
clf;
for i=1:length(beta_list);
    beta = beta_list(i);
    imageplot(phi(f,beta), ['\beta=' num2str(beta)], 2,3,i);
end
