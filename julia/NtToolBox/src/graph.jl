function perform_dijkstra_fm(W, pstart; method="fm", svg_rate=10, niter=Inf, bound="sym")
    """
        # perform_fm_dijstra - *slow* (matlab) implementation of Dijstra and FM
        #
        #   [D,Dsvg,Ssvg] = perform_fm_dijstra(W, pstart, options);
        #
        #   W is an (n,n) metric matrix.
        #   pstart is a (2,k) starting points.
        #   options.method is either 'fm' or 'dijstra'
        #
        #   D is the final distance map to pstart
        #   options.svg_rate gives the rate at wich Dsvg and Ssvg is filled.
        #   options.niter can be used to limit the total number of steps (partial propagation).
        #
        #   Copyright (c) 2012 Gabriel Peyre
    """

##
# Size.

n = size(W,1)

##
# The four displacement vector to go to the four neightbors.

neigh = [[1;0] [-1;0] [0;1] [0;-1]]

##
# For simplicity of implementation, we use periodic boundary conditions.

function symmetrize(x,n)
    if (x<=0)
        x = 2-x;
    elseif (x>n)
        x = 2*n-x
    end
    return x
end

if bound == "per"
    boundary = x -> mod(x-1,n)+1
elseif bound == "sym"
    boundary = x -> [symmetrize(x[1],n), symmetrize(x[2],n)]
else
    error("Works only for per and sym.")
end
#     case 'sym'
#         boundary = @(x)x.*(x<=n & x>0) + (2-x).*(x<=0) + (2*n-x).*(x>n);
#     otherwise
#         error('Works only for per and sym.');
# end

##
# For a given grid index |k|, and a given neighboring index k in \({1,2,3,4}\),
# |Neigh(k,i)| gives the corresponding grid neigboring index.

ind2sub1 = k -> Int.([rem(k-1, n)+1; (k - rem(k-1, n) - 1)/n + 1])
sub2ind1 = u -> Int((u[2]-1)*n + u[1])
Neigh = (k,i) -> sub2ind1(boundary(ind2sub1(k) + neigh[:,i]));
extract   = (x,I) -> x[I]
extract1d = (x,I) -> extract(vec(x),I);

##
# Stack of starting points.
nstart = size(pstart,2)
I = Int.(zeros(nstart))
for i=1:nstart
    I[i] = sub2ind1(pstart[:,i])
end

##
# Initialize the distance to \(+\infty\), excepted for the boundary conditions.

D = zeros(n,n)+Inf; # current distance
for i in 1:nstart
    D[pstart[1,i], pstart[2,i]] = 0;
end

##
# Initialize the state to 0 (unexplored), excepted for the boundary point to \(1\)
# (front).

S = zeros(n,n);
for i in 1:nstart
    S[pstart[1,i],pstart[2,i]] = 1 # open
end

##
# Run!
compt=0
iter = 0
q = 100  # maximum number of saves
Dsvg = zeros(n,n,q)
Ssvg = zeros(n,n,q)
while !(isempty(I)) && iter<=niter
    iter = iter+1
    # pop from stack
    j = sortperm(extract1d(D,I))
    j=j[1]
    i = I[j]
    deleteat!(I,j)
    # declare dead
    u = ind2sub1(i);
    S[u[1],u[2]] = -1;
    # Make a list of neighbors that are not dead
    J = []
    for k in 1:4
        j = Neigh(i,k)
        if extract1d(S,j)!= -1
            # add to the list of point to update
            append!(J,j)
            if extract1d(S,j)==0
                # add to the front
                u = ind2sub1(j)
                S[u[1],u[2]] = 1
                append!(I,j)
            end
        end
    end
    # update neighbor values
    DNeigh = (D,k) -> extract1d(D,Neigh(j,k))
    for j in J
        dx = min(DNeigh(D,1), DNeigh(D,2))
        dy = min(DNeigh(D,3), DNeigh(D,4))
        u = ind2sub1(j)
        w = extract1d(W,j);
        if method=="dijkstra"
            D[u[1],u[2]] = min(dx + w, dy + w)
        else
            Delta = 2*w - (dx-dy).^2
            if (Delta>=0)
                D[u[1],u[2]] = (dx + dy + sqrt(Delta))/ 2
            else
                D[u[1],u[2]] = min(dx + w, dy + w)
            end
        end
    end
    # svd
    t = Base.div(iter,svg_rate)
    if mod(iter,svg_rate)==0 && (t<q)
        Dsvg[:,:,t] = D
        Ssvg[:,:,t] = S
    end
end
return (D,Dsvg,Ssvg)
end
