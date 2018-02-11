gamma = float([x1])
for i in 1:1.5*n/tau
    append!(gamma, gamma[end] - tau*Geval(G,gamma[end]))
    if norm(gamma[end]-x0)<1
        break
    end
end
append!(gamma,x0)
