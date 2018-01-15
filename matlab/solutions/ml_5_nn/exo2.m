Layers = { ...
	[p k],  ... % vanilla logitic regression (convex)
	[p 8 k],  ... % single hidden layer
	[p 3 4 k],  ... % 4 hidden layers
	[p 4 5 8 4 k],  ... % 4 hidden layers
};
tau_list = [.01/10 .01/10 .01/30 .01/40];
clf; 
for il = 1:length(Layers)
    D = Layers{il};
    % layers
    R = length(D)-1;
    A = {}; b = {};
    for r=1:R
        A{r} = randn(D(r+1),D(r));
        b{r} = randn(D(r+1),1);
    end
    % descent
    tau = tau_list(il);
    L = [];
    niter = 12000;
    for it=1:niter
        [L(it),Ag,bg,gx] = EvalNN(x,Y, A,b,loss,rho);
        for r=1:R
            A{r} = A{r} - tau*Ag{r};
            b{r} = b{r} - tau*bg{r};
        end
    end
    % clf; plot(1:niter, L, 'LineWidth', 2); axis tight;
    % probability 
    V = EvalNN(Z,[], A,b,loss,rho);
    U = reshape(SM(V{end}'), [q q k]);
    % same color
    R = zeros(q,q,3);
    for i=1:k
        for a=1:3
            R(:,:,a) = R(:,:,a) + U(:,:,i) .* col(a,i);
        end
    end
    % display
    subplot(2,2,il);
    hold on;
    imagesc(t,t,permute(R, [2 1 3]));
    for j=1:k
        I = find(y==j);
        plot(x(1,I), x(2,I), '.', 'color', col(:,j)*.8, 'MarkerSize', 20);
    end
    axis equal; axis off;
    title(['D=' num2str(D)]);
end
%EXO