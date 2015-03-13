---
layout: page
title: "Installation"
description: "for Julia"
header-img: "img/hokusai-5.jpg"
---

Installation of the toolboxes packages
------------

Each tour makes use some of a toolbox package, that needs to be [downloaded as a .zip file](https://github.com/gpeyre/numerical-tours/raw/master/julia/nt_toolbox.zip), and then unzipped within your working directory (from which you will run the notebooks): Alternatively, you can download the whole [numerical_tours][1].

You will also need:

* [PyPlot](https://github.com/stevengj/PyPlot.jl) package to be able to call [Matplotlib](http://matplotlib.org/) for displaying plots.
* [Autoreload](https://github.com/malmaud/Autoreload.jl) package to ease developpement (automatic reloading of modified external files).

You can install these package using the

> Pkg.add("PyPlot")

and

> Pkg.add("Autoreload")

from the Julia command line.

Installation of IPython
------------

If you intend to run the tours are IPython notebook (which we recommend), you need to install [IPython notebook][2] to run the code.
You will also need to install [IJulia](https://github.com/JuliaLang/IJulia.jl).


Using IPython
------------

To run the Julia's tours as IPython notebook using IJulia, you need to call from a terminal (from the correct directory location):

> ipython notebook --profile julia


[1]: https://github.com/gpeyre/numerical-tours/archive/master.zip
[2]: http://ipython.org/install.html
[3]: http://arokem.github.io/python-matlab-bridge/
