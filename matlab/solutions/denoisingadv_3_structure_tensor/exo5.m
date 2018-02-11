[X,Y] = meshgrid(1:n,1:n);
[XI,YI] = meshgrid(1:1/2:n+1/2,1:1/2:n+1/2); XI(:,end) = n; YI(end,:) = n;
f1 = interp2(X,Y,f,XI,YI);
% Compute the gradient field.
G = grad(f1);
% Compute the rank-1 tensor field associated to the gradient rotated by
% 90°.
T0 = zeros(n*2,n*2,2,2);
T0(:,:,1,1) = G(:,:,2).^2;
T0(:,:,2,2) = G(:,:,1).^2;
T0(:,:,1,2) = -G(:,:,1).*G(:,:,2);
T0(:,:,2,1) = -G(:,:,1).*G(:,:,2);
% Blur the tensor field.
sigma = 20;
T = perform_blurring(T0, sigma);
% Sub-sample it.
T = T(1:2:end,1:2:end,:,:);
options.sub = 4;
clf;
plot_tensor_field(T, [], options);
