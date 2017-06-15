function perform_haar_transf(M, Jmin, dir)
    n = size(M)[1]
    Mf = copy(M)
    if dir == 1
        MW = copy(M)

        for j in Jmin : Int(log2(n))
            p = Int(n/2^(j-1))
            sel = 1 : p
            even = 1 : 2 : p
            odd = 2 : 2 : p
            # average/ difference along X
            MW[sel, sel, sel] = cat(1, (MW[even, sel, sel] + MW[odd, sel, sel])./sqrt(2), (MW[even, sel, sel] - MW[odd, sel, sel])./sqrt(2))
            # average/ difference along Y
            MW[sel, sel, sel] = cat(2, (MW[sel, even, sel] + MW[sel, odd, sel])./sqrt(2), (MW[sel, even, sel] - MW[sel, odd, sel])./sqrt(2))
            # average/ difference along Z
            MW[sel, sel, sel] = cat(3, (MW[sel, sel, even] + MW[sel, sel, odd])./sqrt(2), (MW[sel, sel, even] - MW[sel, sel, odd])./sqrt(2))
            Mf = copy(MW)
        end
    else

        M1 = copy(M)

        for j in Int(log2(n)) : -1 : Jmin
            p = Int(n/2^j)
            sel = 1 : p
            sel1 = 1:2*p
            selw = p+1:2*p
            even = 1:2:2*p
            odd = 2:2:2*p
            # average/ difference along X
            A = M1[sel, sel1, sel1]
            D = M1[selw, sel1, sel1]
            M1[even, sel1, sel1] = (A + D)./sqrt(2)
            M1[odd, sel1, sel1] = (A - D)./sqrt(2)
            # average/ difference along Y
            A = M1[sel1, sel, sel1]
            D = M1[sel1, selw, sel1]
            M1[sel1, even, sel1] = (A + D)./sqrt(2)
            M1[sel1, odd, sel1] = (A - D)./sqrt(2)
            # average/ difference along Z
            A = M1[sel1, sel1, sel]
            D = M1[sel1, sel1, selw]
            M1[sel1, sel1, even] = (A + D)./sqrt(2)
            M1[sel1, sel1, odd] = (A - D)./sqrt(2)
            Mf = copy(M1)
        end
    end

    return Mf
end
