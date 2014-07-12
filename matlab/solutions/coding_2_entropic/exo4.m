e = -sum( h.*log2(h) );
% compute the differencial of coding for a varying length signal
err = []; 
slist = 4:12;
for i = 1:length(slist)
    n = 2^slist(i);
    x = rand_discr(h, n);
    % coding
    y = perform_arith_fixed(x(:),h);
    nb = length(y);
    e1 = nb/n; % number of bits per symbol
    err(i) = e1 - e;
end
clf;
plot(slist, err, '.-'); axis('tight');
set_label('log2(size)', '|entropy-nbr.bits|');
