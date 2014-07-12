pathsD = {};
for k=1:nend
    % path purely on edges
    vlist = pend(k);
    vprev = D(vlist(end));
    while true
        x0 = vlist(end);
        r = vring{x0};
        [v,J] = min(D(r));
        x = r(J);
        if v>=vprev || v==0
            break;
        end
        vprev = v;
        vlist(end+1) = x;
    end
    pathsD{end+1} = vertex(:,vlist);
end
plot_fast_marching_mesh(vertex,faces, Q, pathsD, options); 
