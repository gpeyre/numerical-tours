s1 = [2:nbins; nbins-1]
s2 = [2; 1:nbins-1]
I = find( (h[s1].<h) & (h[s2].<h) )
# sort in descending order
v = sort(h[I])
u = sortperm(h[I])
if v[1]<v[end]
    u = reverse(u)
end
theta1 = t[I[u[1:3]]]
theta1 = theta1[:]'
M1 = [cos(theta1); sin(theta1)]
println("--- M ---")
println(M)
println("--- M1 ---")
println(M1)
