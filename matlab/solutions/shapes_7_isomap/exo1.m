for i=1:n
    % progressbar(i,n);
    D = min(D,repmat(D(:,i),[1 n])+repmat(D(i,:),[n 1])); 
end
