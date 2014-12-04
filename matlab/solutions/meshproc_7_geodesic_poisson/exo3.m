nlist = [1 2 3 4]; 
nlist = [3 6 10 20]; 
k = 1;
i = 21000; % for elephant (nose)
I = [i];
delta = zeros(n,1); 
delta(i) = 1;
t = 1/1000; 
clf;
for i=1:max(nlist)
    % solve
    u = (Ac+t*DeltaCot)\delta;
   	% re-integrate normalized gradient
    phi = -Delta \ Div( normalize( Grad(u) ) );
    % display
    if k<=4 && i==nlist(k)
        subplot(2,2,k); hold on;
        plot3(X(1,I), X(2,I), X(3,I), 'r.', 'MarkerSize', 25);
        options.face_vertex_color = DispFunc(phi);
        plot_mesh(X,F,options);
        drawnow;
        axis('tight');
        colormap parula(256);
        k = k+1;
    end
    % add fatherest point
    [~,j] = max(phi);
    delta(j) = 1;
    I(end+1) = j;
end
