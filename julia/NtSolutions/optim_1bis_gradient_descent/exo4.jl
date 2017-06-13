function ProjGradDescent(Grad_f, Proj_Omega, y, mask, x0, nbiter, tau)
    x = x0;
    for iter in 1:nbiter  
        x = Proj_Omega(x - tau*Grad_f(x), y, mask)  
    end
    return x
end;