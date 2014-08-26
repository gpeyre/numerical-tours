slist = linspace(1e-9, 1, 1000);
RL = [];
for is=1:length(slist)
    % observations
    sigma = slist(is);
    y = y0 + sigma*norm(y0)*w;
    % SVD
    [U,S,V] = svd(MusicHankel(y),0); S = diag(S);
    Ubot = U(:,N+1:end);
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
    % position
    [x1,I] = sort(mod(angle(R),2*pi)/(2*pi));
    R = R(I);
    % Keep only the best N ones.
    if 0
        [~,I] = sort(abs(abs(R)-1));
        x1 = x1(I(1:N));
    end
    RL(:,end+1) = R;
end
%
clf; hold on;
plot(exp(2i*pi*z), 'k');
%
cm = jet(length(slist));
for is=1:length(slist)
    plot(RL(:,is), '.',  'MarkerSize', ms, 'color', cm(is,:));
end
%
axis equal; box on;
axis([-1 1 -1 1]*1.1);
axis off;
