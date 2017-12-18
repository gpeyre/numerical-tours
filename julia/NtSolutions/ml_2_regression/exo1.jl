q = 50
lambda_list = linspace(0,20,q)
W = reshape([], size(X0,2),0) # empty array
E = []
for i=1:q
    lambda = lambda_list[i]
    w = (X0'*X0+lambda*eye(p)) \ (X0'*y0)
    W = hcat(W,w) # bookkeeping
    append!(E,sqrt( sum( (X1*w-y1).^2 ) / n1 ))
end
# Display error evolution.
plot(lambda_list, E, lw=2)
axis("tight")
xlabel("\lambda")
ylabel("E");
