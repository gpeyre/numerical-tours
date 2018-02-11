% compute the energy of each point
d = sum(P.^2,2);
% extract only a small sub-set
rho = .1;
[v,I] = sort(d);
if v(1)<v(length(v))
    I = reverse(I);
end
I = I(1:round(rho*length(I)));
P1 = P(I,:);  % transformed points
% compute histogram
Theta = mod(atan2(P1(:,2),P1(:,1)), pi());
% display histograms
nbins = 200;
[h,t] = hist(Theta, nbins);
h = h/sum(h);
clf;
bar(t,h); axis('tight');
