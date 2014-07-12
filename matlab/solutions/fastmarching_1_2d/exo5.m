% Hint: you can use a convolution |U=conv2(double(Q==2),h,'same')| with a
% well chose kernel |h| to located the points |U>0| with at least 1
% neighbor.
h = [0 1 0; 1 0 1; 0 1 0];
B = (Q==1) & (conv2(double(Q==2),h,'same')>0);
U = find(B);
[xa,xb] = ind2sub(size(f),U);
[xa,I] = sort(xa); xb = xb(I); U = U(I);
dU = D(U);
k = [65 259];
x1 =[ [xa(k(1));xb(k(1))] [xa(k(2));xb(k(2))] ];
%
clf;
% subplot(2,1,1);
hold on;
imageplot(A, 'Boundary points');
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);
h = plot(xb,xa, 'g'); set(h, 'LineWidth', 2);
h = plot(x1(2,:), x1(1,:), '.b'); set(h, 'MarkerSize', 25);
%
if 0
subplot(2,1,2);
hold on;
h = plot( [k(1) k(1)], [min(dU) max(dU)], 'r:' ); set(h, 'LineWidth', 2);
h = plot( [k(2) k(2)], [min(dU) max(dU)], 'r:' ); set(h, 'LineWidth', 2);
h = plot(dU); axis('tight');title('D along the boundary'); set(h, 'LineWidth', 2);
end
