function perform_thresholding(f, M, t)
    """
        Only 3 types of thresholding currently implemented
    """
    if t == "largest"
        a = sort(abs(f)[:])[end:-1:1] #sort a 1D copy of F in descending order
        T = a[M]
        y = f.*(abs(f) .> T)
    elseif t == "soft"
        s = abs(f) - M
        s = (s + abs(s))./2
        y = (f./abs(f)).*s
    elseif t == "hard"
        y = f.*(abs(f) .> M)
    end
    return y
end
