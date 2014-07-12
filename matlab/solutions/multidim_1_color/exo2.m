a = @(g)g(:,:,2) .* cos(g(:,:,1));
b = @(g)g(:,:,2) .* sin(g(:,:,1));
% This ugly code is for Scilab compatibility
c = @(g)cat(3, g(:,:,3), a(g), b(g));
hsv12rgb = @(g)applymat(c(g),T);
clf;
theta = linspace(0,pi/2,6);
for i=1:length(theta)
    g1 = g;  g1(:,:,1) = g1(:,:,1) + theta(i); 
    imageplot(clamp(hsv12rgb(g1)), ['\theta=' num2str(theta(i))], 2,3,i);
end
