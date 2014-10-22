slist = round([P/16 P/8 P/4 P/2]);
slist = round([2 3 6 10]);
clf;
for it_s=1:length(slist)
    s = slist(it_s);
    %
    I = randperm(N); I = I(1:s);
    x0 = zeros(N,1); 
    % x0(1:s) = 1;
    x0(I) = sign(randn(s,1));
    sigma = 0.1 * norm(Phi*x0)/sqrt(P);
    y = Phi*x0 + sigma*randn(P,1);
    %
    [X,Lambda] = perform_homotopy(Phi,y);
    % Fuchs criteria
    J = setdiff(1:N, I);
    U = Phi(:,J)' * pinv(Phi(:,I))';
    c = U * sign(x0(I));
    Fuchs = max(abs(c));
    % ERC 
    ERC = norm(U, 'inf');    
    %
    subplot(2,2,it_s);
    title(['s=' num2str(s), ', F=' num2str(Fuchs, 2) ', ERC=' num2str(ERC,2)]);
    I = find(x0~=0); J = setdiff(1:N,I);
 	hold on;
    h = plot(Lambda, X(I,:)', '-'); set(h, 'LineWidth', 2);
    h = plot(Lambda, X(J,:)', 'k-'); set(h, 'LineWidth', 2);
    axis tight;
    box on;
end
