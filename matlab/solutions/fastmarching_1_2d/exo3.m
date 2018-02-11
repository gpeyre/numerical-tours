elist = [.1 .5 1 10];
clf;
for ie=1:length(elist)
    epsilon = elist(ie);
    W = epsilon + abs(f-f(x0(1),x0(2)));
    [D,S] = perform_fast_marching(1./W, x0, options);
    G0 = grad(D, options);
    G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);
    % 
    gamma = x1;
    for i=1:1.5*n/tau
        gamma(:,end+1) = gamma(:,end) - tau*Geval(G, gamma(:,end));
        if norm(gamma(:,end)-x0)<1
            break;
        end
    end
    gamma(:,end+1) = x0;
    %
    subplot(2,2,ie); hold on;
    imageplot(f);
    h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
    h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
    h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij; 
    title(['\epsilon=' num2str(epsilon)]);
end
