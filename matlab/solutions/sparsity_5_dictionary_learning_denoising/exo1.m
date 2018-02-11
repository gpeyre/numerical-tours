tx = zeros(p,size(y,2));
tu = D*tx;
E = [];
E1 = [];
for i=1:niter
    tx1 = tx; tu1 = tu;
    if 1
        %%%% MODE 1 %%%%
        % rProxG
        [tx2,tu2] = ProxG(tx1, tu1, gamma);
        tx1 = 2*tx2 - tx1; tu1 = 2*tu2 - tu1;
        % rProxF
        [tx2,tu2] = ProxF(tx1, tu1, gamma);
        tx1 = 2*tx2 - tx1; tu1 = 2*tu2 - tu1;
        % Average
        tx = (1-mu/2)*tx + mu/2*tx1;
        tu = (1-mu/2)*tu + mu/2*tu1;
        % Converging sequence
        [x,u] = ProxG(tx,tu,gamma);
    else
        %%%% MODE 2 %%%%
        % rProxF
        [tx2,tu2] = ProxF(tx1, tu1, gamma);
        tx1 = 2*tx2 - tx1; tu1 = 2*tu2 - tu1;
        % rProxG
        [tx2,tu2] = ProxG(tx1, tu1, gamma);
        tx1 = 2*tx2 - tx1; tu1 = 2*tu2 - tu1;
        % Average
        tx = (1-mu/2)*tx + mu/2*tx1;
        tu = (1-mu/2)*tu + mu/2*tu1;
        % Converging sequence
        [x,u] = ProxF(tx,tu,gamma);
    end
    % 
    E(end+1,:) = sum(abs(x));
    E1(end+1,:) = sqrt(sum( (D*x-y).^2, 1 ))/epsilon;
end
[~,I] = sort(sqrt(sum(Y.^2)), 'descend');
sel = I(1:5);
clf;
subplot(2,1,1);
plot(log10(E(1:end/2,sel) - repmat(min(E(:,sel),[],1),[niter/2 1])  ));
axis tight;
title('log_{10}(|x|_1 - |x^*|_1)');
subplot(2,1,2);
plot(log(abs(E1(:,sel)-1)));
axis tight;
title('log_{10}| |Dx-y| - \epsilon|');
