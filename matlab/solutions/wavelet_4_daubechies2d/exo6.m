fWSep = f;
Jmin = 1;
for i=1:n
    fWSep(:,i) = perform_wavortho_transf(fWSep(:,i),Jmin,+1);
end
for i=1:n
    fWSep(i,:) = perform_wavortho_transf(fWSep(i,:)',Jmin,+1)';
end
