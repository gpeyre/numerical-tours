% entropy bound
e = -sum( log2(h).*h );
disp(['Entropy=' num2str(e) '.']);
qlist = 1:10;
err = [];
for q=qlist
    % lifting
    n1 = ceil(n/q)*q;
    x1 = x;
    x1(length(x1)+1:n1) = 1;
    x1 = reshape(x1,[q n1/q]);
    [Y,X] = meshgrid(1:n1/q,0:q-1);
    x1 = sum( (x1-1) .* (m.^X), 1 )' + 1;
    % Probability table
    H = h;
    for i=1:q-1
        H = kron(H,h);
    end
    % compute the tree
    T = compute_hufftree(H);
    % do the coding
    y = perform_huffcoding(x1,T,+1);
    % average number of bits
    e1 = length(y)/length(x);
    err(q) = e1-e;
    disp(['Huffman(block size ' num2str(q) ')=' num2str(e1)]);
end
clf;
plot(qlist,err, '.-');
set_label('q', 'entropy-code.length');
axis('tight');
