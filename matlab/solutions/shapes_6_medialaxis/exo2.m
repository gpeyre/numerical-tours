thresh = [10 20 50 100];
clf;
for i=1:4
    subplot(2,2,i);
    B = G>thresh(i);
    A = M;
    A(B==1) = 0;
    imageplot(-A);
end
