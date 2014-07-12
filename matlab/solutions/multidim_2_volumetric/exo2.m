M1 = MWT;
for j=log2(n):-1:1
    p = n/2^(j);
    sel = 1:p; 
    sel1 = 1:2*p;
    selw = p+1:2*p;
    % average/difference along X
    A = M1(sel,  sel1, sel1);
    D = M1(selw, sel1, sel1);
    M1(1:2:2*p,sel1,sel1) = (A+D)/sqrt(2);
    M1(2:2:2*p,sel1,sel1) = (A-D)/sqrt(2);
    % average/difference along Y
    A = M1(sel1, sel,  sel1);
    D = M1(sel1, selw, sel1);    
    M1(sel1,1:2:2*p,sel1) = (A+D)/sqrt(2);
    M1(sel1,2:2:2*p,sel1) = (A-D)/sqrt(2);
    % average/difference along Z
    A = M1(sel1, sel1, sel);
    D = M1(sel1, sel1, selw);    
    M1(sel1,sel1,1:2:2*p) = (A+D)/sqrt(2);
    M1(sel1,sel1,2:2:2*p) = (A-D)/sqrt(2);
end
