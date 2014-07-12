P = 0;
for i=1:length(S)
    P = P + length(S{i}(:));
end
disp( strcat(['True redundancy of the dictionary=' num2str(P/n) '.']) );
