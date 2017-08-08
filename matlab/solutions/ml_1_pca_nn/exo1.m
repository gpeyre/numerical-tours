D = distmat(X0,X1);
[Ds,I] = sort(D,1);
Y = y(I);
Rmax = 50;
for R=1:Rmax
    if R==1
        C = Y(1,:); 
    else
        h = hist(Y(1:R,:), 1:k);
        [~,C] = max(h);
    end
    % correct classification 
    S(R) = sum(C(:)==y1)/p1;
end
clf;
plot(1:Rmax, S, '.-', 'MarkerSize', ms);
axis tight;
SetAR(1/2);
xlabel('R'); ylabel('S');
