J = log2(N)-1;
Jmin = 0;
fw = f;
clf;
subplot(4,1,1);
plot(f); axis('tight'); title('Signal');
for j=J:-1:Jmin
    fw(1:2^(j+1)) = haar(fw(1:2^(j+1)));

    
    j1 = J-j;
    if j1<3
        d = fw(2^j+1:2^(j+1));
        subplot(4,1,j1+2);
        plot(1:2^(j1+1):N,d);  axis('tight');
        title( strcat(['Details, j=' num2str(j)]) );
    end   
    
    if 0
    A = fw(1:2^(j+1));
    Coarse = ( A(1:2:length(A)) + A(2:2:length(A)) )/sqrt(2);
    Detail = ( A(1:2:length(A)) - A(2:2:length(A)) )/sqrt(2);
    A = cat3(1, Coarse, Detail );
    fw(1:2^(j+1)) = A;
    j1 = J-j;
    if j1<3
        subplot(4,1,j1+2);
        plot(1:2^(j1+1):N,Detail);  axis('tight');
        title( strcat(['Details, j=' num2str(j)]) );
    end   
    end
end
