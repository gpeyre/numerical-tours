def exo1():
    """
    Try with other alphamapping and colormapping
    et up a colormap
    et up an alpha map
    efresh the rendering
    """
    h = vol3d('cdata', M, 'texture', '2D')
    view(3); axis off
    colormap bone(256)
    options.sigma = .08; % control the width of the non-transparent region
    options.center = .6; % here a value in [0, 1]
    a = compute_alpha_map('gaussian', options); % you can plot(a) to see the alphamap
    vol3d(h)


def exo2():
    """
    select the point (x,y) of minimum value in the slice |D(:,:,n-delta)|.
    hint: use functions 'min' and 'ind2sub'
    """
    d = D(: , : , n-delta)
    [tmp, I] = min(d(: ))
    [x, y] = ind2sub([n n], I(1))
    end_point = [x; y; n-delta]


def exo3():
    """
    Select other starting points. In order to do so, ask the user to
    click on a starting point in a given horizontal slice |W(:,:,delta)|.
    You can do this by using |ginput(1)|
    on the plane |Z=delta|.
    """
    if 0
        clf; imageplot(M(: , : , delta))
        title('Pick starting point')
        start_point = round(ginput(1))
        start_point = [start_point(2); start_point(1); delta]


