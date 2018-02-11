gamma = x1;
for i=1:1.5*n/tau
    gamma(:,end+1) = gamma(:,end) - tau*Geval(G, gamma(:,end));
    if norm(gamma(:,end)-x0)<1
        break;
    end
end
gamma(:,end+1) = x0;
