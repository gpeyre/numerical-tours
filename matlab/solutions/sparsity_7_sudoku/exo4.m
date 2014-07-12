Ablock = [];
for i=1:p
    for j=1:p
        for k=1:n          
            Z = zeros(n,n,n);
            Z((1:p) + (i-1)*p,(1:p) + (j-1)*p,k) = 1;
            Ablock(end+1,:) = Z(:)';
        end        
    end
end
clf;
imageplot(Ablock);
