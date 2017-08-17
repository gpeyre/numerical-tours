function rand_discr(p; m = 1)
    """
        rand_discr - discrete random generator

          y = rand_discr(p, n);

          y is a random vector of length n drawn from
          a variable X such that
              p(i) = Prob( X=i )

          Copyright (c) 2004 Gabriel Peyré
    """

    # makes sure it sums to 1
    p = p/sum(p)

    n = length(p)
    coin = rand(m)
    cumprob = [0.]
    append!(cumprob, cumsum(p))
    sample = zeros(m)

    for j in 1:n
        ind = find((coin .> cumprob[j]) & (coin .<= cumprob[j+1]))
        sample[ind] = j
    end
    return sample
end
