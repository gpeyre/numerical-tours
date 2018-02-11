% for bookeeping
X = []; Lambda = []; Sparsity = [];
% initialization
C = Phi'*y;
[lambda,I] = max(abs(C));
x = zeros(N,1);
X(:,end+1) = x; Lambda(end+1) = lambda; Sparsity(end+1) = length(I);
clf; 
niter = 1000;
for it=1:niter
    % complementary support
    J = setdiff(1:N, I);
    % correlations
    c = Phi'*(y-Phi*x);
    if norm(c)<1e-6
        break;
    end
    if it<=6
        % Display to check that x_lambda is indeed a solution
        subplot(2,3,it);
        hold on;
        stem( I, c(I)/lambda, 'b.' );
        stem( J, c(J)/lambda, 'r.' );
        plot([1 N], [1 1], 'k--');
        plot([1 N],-[1 1], 'k--');
        axis([1 N -1.05 1.05]);
        title(['Step ' num2str(it)]);
    end
    % update direction
    d = zeros(N,1);
    d(I) = (Phi(:,I)'*Phi(:,I)) \ sign(c(I));
    % useful vector
    v = Phi(:,I)*d(I);
    % Compute minimum |gamma| so that situation 1) is in force.
    w = ( lambda-c(J) ) ./ ( 1 - Phi(:,J)'*v );
    gamma1 = min(w(w>0)); 
    if not(isempty(gamma1))
        i1 = J( w==gamma1 );
    end
    % Compute minimum |gamma| so that situation 2) is in force.
    w = ( lambda+c(J) ) ./ ( 1 + Phi(:,J)'*v );
    gamma2 = min(w(w>0));  
    if not(isempty(gamma2))
        i2 = J( w==gamma2 );
    end
    % Compute minimum |gamma| so that situation 3) is in force.
    w = -x(I)./d(I);
    gamma3 = min(w(w>0));     
    if not(isempty(gamma3))
        i3 = I( w==gamma3 );
    end
    % any condition is in force
    gamma = min([gamma1 gamma2 gamma3]);
    if isempty(gamma)
        break;  % final solution reached
    end
    % new solution
    x = x + gamma*d;
    lambda = lambda - gamma;
    if gamma==gamma1
        I = [I i1]; 
    elseif gamma==gamma2
        I = [I i2]; 
    elseif gamma==gamma3
        I(I==i3) = [];
        x(i3) = 0;
    end
    % record sparsity and lambda
    X(:,end+1) = x; Lambda(end+1) = lambda;
    Sparsity(end+1) = length(I);
end
