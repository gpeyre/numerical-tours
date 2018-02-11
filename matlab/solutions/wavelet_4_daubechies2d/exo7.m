f1 = fWSep;
for i=1:n
    f1(:,i) = perform_wavortho_transf(f1(:,i),Jmin,-1);
end
for i=1:n
    f1(i,:) = perform_wavortho_transf(f1(i,:)',Jmin,-1)';
end
