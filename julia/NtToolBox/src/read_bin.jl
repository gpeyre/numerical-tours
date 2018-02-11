function read_bin(name, ndims = 2)
    """
        reading from a binary file in 3 dimensions
    """

    file = open(name)
    n = Int(read(file, UInt16, 1)[1])
    p = Int(read(file, UInt16, 1)[1])
    q = 1
    if ndims == 3
        q = Int(read(file, UInt16, 1)[1])
    end
    M = read(file, UInt8, n*p*q)
    M = reshape(M, (n, p, q))
    close(file)
    return M
end
