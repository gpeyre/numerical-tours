% control mesh.
f0 = rand(12,3);
Jmax = 4; ms = 20; lw = 1.5;
f = f0;
for j=0:Jmax
    f = cat(2, upsampling(f(:,1)), upsampling(f(:,2)), upsampling(f(:,3))  );
    f = cat(2, cconvol(f(:,1),h), cconvol(f(:,2),h), cconvol(f(:,3),h) );    
end
clf;
subplot(1,2,1);
hold on;
hh = plot3([f(:,1);f(1,1)], [f(:,2);f(1,2)], [f(:,3);f(1,3)], 'k-');
set(hh, 'MarkerSize', ms);
set(hh, 'LineWidth', lw); 
hh = plot3([f0(:,1);f0(1,1)], [f0(:,2);f0(1,2)], [f0(:,3);f0(1,3)], 'r.--');
set(hh, 'LineWidth', lw);
axis('tight'); box('on'); view(3);  % axis('off');
subplot(1,2,2);
hold on;
hh = plot3([f(:,1);f(1,1)], [f(:,2);f(1,2)], [f(:,3);f(1,3)], 'k-');
set(hh, 'MarkerSize', ms);
set(hh, 'LineWidth', lw); 
hh = plot3([f0(:,1);f0(1,1)], [f0(:,2);f0(1,2)], [f0(:,3);f0(1,3)], 'r.--');
set(hh, 'LineWidth', lw);
axis('tight'); box('on');
view(70,25);
