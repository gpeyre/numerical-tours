% Coefficients
B = [];
for j=1:size(Ubot,2)
    u = Ubot(:,j);
    v = flipud(conj(u)); 
    B(:,j) = conv(u,v);
end
C = sum(B,2);
% Roots
R = roots(C(end:-1:1));
% keep those inside
R = R(abs(R)<=1);
% display
clf; hold on;
plot(exp(2i*pi*z), 'k');
plot(R, 'b.', 'MarkerSize', ms);
axis equal; box on;
axis([-1 1 -1 1]*1.1);
axis off;
