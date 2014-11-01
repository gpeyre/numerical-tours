f1 = fW;
clf;
for j=Jmin:Jmax
    A = f1(1:2^(j+1),1:2^(j+1));
    for d=1:2
        if d==1
            Coarse = A(1:2^j,:);
            Detail = A(2^j+1:2^(j+1),:);
        else
            Coarse = A(:,1:2^j);
            Detail = A(:,2^j+1:2^(j+1));                
        end
        Coarse = cconvol(upsampling(Coarse,d),reverse(h),d);
        Detail = cconvol(upsampling(Detail,d),reverse(g),d);
        A = Coarse + Detail;
        j1 = Jmax-j;
        if j1>0 & j1<5
            subplot(2,2,j1);
            imageplot(A, strcat(['Partial reconstruction, j=' num2str(j)]));
        end
    end
    f1(1:2^(j+1),1:2^(j+1)) = A;
end
