% test for mesh subdivision

path(path, '../toolbox_graph_data/off/');
path(path, '../toolbox_graph_data/off/low_poly/');

name = 'ico';
name = 'mushroom';
name = 'venus';
name = 'mannequin';
name = 'patch';

vertex = {}; face = {};
[vertex{1},face{1}] = read_mesh(name);
options.name = name;

J = 3; %  number of scales
sublist = {'sqrt3', 'loop', 'butterfly', 'linear4'};
sublist = {'sqrt3'};

rep = 'results/subdivision-surfaces/';
if not(exist(rep))
    mkdir(rep);
end

display_control = 0;
for it=1:length(sublist)
    options.sub_type = sublist{it};
    disp(['--> Testing ' sublist{it} ' subdivision.']);
    for j=1:J
        if j>1
            [vertex{j},face{j}] = perform_mesh_subdivision(vertex{j-1}, face{j-1}, 1, options);
        end
        clf;
        hold on;
        plot_mesh(vertex{j},face{j},options);
        shading faceted;
        if display_control
            h = plot3(vertex{1}(1,:), vertex{1}(2,:), vertex{1}(3,:), 'r.'); % control mesh
            set(h, 'MarkerSize', 20);
            hold off;
        end
        saveas(gcf, [rep name '-subdivision-' sublist{it} '-' num2str(j) '.png'], 'png');
    end
end