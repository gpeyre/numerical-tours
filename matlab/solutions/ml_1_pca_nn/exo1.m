D = distmat(X0,X1);
[Ds,I] = sort(D,1);
ys = y(I);
Rmax = 50;
S = [];
for R=1:Rmax
    if R==1
        C = ys(1,:); 
    else
        h = hist(ys(1:R,:), 1:k);
        [~,C] = max(h);
    end
    % correct classification 
    S(R) = sum(C(:)==y1)/n1;
end
clf;
% plot(1:Rmax, S, '.-', 'MarkerSize', ms);
bar(1:Rmax, S);
axis tight;
axis([1 Rmax min(S)*.99 1]);
SetAR(1/2);
xlabel('R'); ylabel('S');
