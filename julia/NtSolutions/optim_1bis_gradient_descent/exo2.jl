function GradDescent(Grad_f, y, Lambda, x0, nbiter, tau)
    # put your code here.
    x = x0;
    for iter in 1:nbiter  
        x = x - tau * Grad_f(x, y, Lambda)  
    end
    return x
end;