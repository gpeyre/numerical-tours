tau = .01/10;
if R==5
    tau = .01/50;
end
if R==3
    tau = .01/80;
end
L = [];
niter = 8000;
for it=1:niter
    [L(it),Ag,bg,gx] = EvalNN(x,Y, A,b,loss,rho); 
    for r=1:R
        A{r} = A{r} - tau*Ag{r};
        b{r} = b{r} - tau*bg{r};
    end
end
clf;
plot(1:niter, L, 'LineWidth', 2);
xlabel('iter');  ylabel('L');
axis tight;
set(gca, 'FontSize', 15);
