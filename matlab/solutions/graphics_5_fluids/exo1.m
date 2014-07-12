rho_list = [2 4 8 16];
clf;
for i=1:length(rho_list)
    rho = rho_list(i);
    imageplot(W(f,rho*U), ['\rho=' num2str(rho)], 2,2,i);
end
