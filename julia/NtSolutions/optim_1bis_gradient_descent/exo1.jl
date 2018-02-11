function Grad_f(x, y, Lambda)
    # put your code here.
    return x - y + Lambda * Dadj(D(x))
end;