def exo1():
    """
    Histogram equalization is an orthogonal projector that maps the values
    of one signal onto the values of the other signal. This is achieved by
    assiging the sorted of ont signal to the sorted values of the other
    signla. Implement this for the two images.
    """
    [v1, I1] = sort(M1(: ))
    [v2, I2] = sort(M2(: ))
    Meq = M1
    Meq(I1) = v2
    Meq = reshape(Meq, [n n])
    imageplot(M1, 'Original', 1, 2, 1)
    imageplot(Meq, 'Equalized', 1, 2, 2)


