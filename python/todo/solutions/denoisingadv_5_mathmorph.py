def exo1():
    """
    Display structuring elements of increasing sizes.
    """
    clf; i = 1
    for w in [1 1.5 3 7]:
        imageplot(strel(w), ['w = ' num2str(w)], 1, 4, i)
        i = i + 1


def exo2():
    """
    Test with structing elements of increasing size.
    """
    i = 0
    for w  in  [1 1.5 2 4]:
        i = i + 1
        imageplot(dillation(M, w), ['w = ' num2str(w)], 2, 2, i)


def exo3():
    """
    Test with structing elements of increasing size.
    """
    i = 0
    for w  in  [1 1.5 2 4]:
        i = i + 1
        imageplot(errosion(M, w), ['w = ' num2str(w)], 2, 2, i)


def exo4():
    """
    Test with structing elements of increasing size.
    """
    i = 0
    for w  in  [1 1.5 2 4]:
        i = i + 1
        imageplot(opening(M, w), ['w = ' num2str(w)], 2, 2, i)


