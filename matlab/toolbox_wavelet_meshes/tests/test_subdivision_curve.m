% test for subdivision curves


options.h = [1 4 6 4 1];

rep = 'results/subdivision-curve/';
if not(exist(rep))
    mkdir(rep);
end


%%%%%%% subdivision function %%%%%%%%
name = 'rand';
name = 'dirac';
name = 'curve';
name = 'square';

switch name
    case 'rand'
        n0 = 9;
        f0 = rescale(rand(1,n0), .05,.95);
    case 'dirac'
        f0 = [0 0 0 1 0 0];
    case 'curve'
        f0 = []; b = 1;
        while b==1
            clf;
            if size(f0,2)>1
                plot(f0(1,:), f0(2,:), '.-');
            end
            axis([0 1 0 1]); box on;
            [x,y,b] = ginput(1);
            f0(:,end+1) = [x;y];
        end
    case 'square'
        f0 = [0 0 1 1; 0 1 1 0];
        f0 = rescale(f0,.05,.95);
end
n0 = size(f0,2);
x0 = linspace(0,1,n0+1);

Jmax = 5; ms = 20; lw = 1.5;
for j=1:Jmax
    f = perform_curve_subdivision(f0, j, options);
    x = linspace(0,1,size(f,2)+1);
    clf;
    hold on;
    if size(f0,1)>1
        h = plot([f(1,:) f(1,1)], [f(2,:) f(2,1)], 'k.-');
    else
        h = plot(x, [f f(1)], 'k.-');
    end
    set(h, 'MarkerSize', ms);
    set(h, 'LineWidth', lw);
    if size(f0,1)>1
        h = plot([f0(1,:) f0(1,1)],[f0(2,:) f0(2,1)], 'r.--');
    else
        h = plot(x0,[f0 f0(1)], 'r.--');
    end
    set(h, 'LineWidth', lw);
    hold off;
    axis([0 1 0 1]); box on;
    saveas(gcf, [rep 'subdivision-func-' name '-' num2str(j) '.png'], 'png');
end
