---
layout: page
title: "Installation"
description: "for Matlab"
header-img: "img/hokusai-5.jpg"
---

Installation of the toolboxes
------------

Each tour makes use some of the following toolboxes, that needs to be downloaded as .zip file, and then unzipped within your working directory (from which you will run the notebooks):

* [toolbox_general](https://github.com/gpeyre/numerical-tours/raw/master/matlab/toolbox_general.zip)
* [toolbox_signal](https://github.com/gpeyre/numerical-tours/raw/master/matlab/toolbox_signal.zip)
* [toolbox_graph](https://github.com/gpeyre/numerical-tours/raw/master/matlab/toolbox_graph.zip) (needed only for a few tours)
* [toolbox_wavelet_meshes](https://github.com/gpeyre/numerical-tours/raw/master/matlab/toolbox_wavelet_meshes.zip) (needed only for a few tours)

Alternatively, you can download the whole [numerical_tours][1].


Scilab and Octave
------------

Note that a lot of Matlab tours are also compatible with [Scilab](http://www.scilab.org/) and with [GNU Octave](https://www.gnu.org/software/octave/).


Installation of IPython
------------

If you intend to run the tours are IPython notebook (which we recommend), you need to install IPython [notebook][2] to run the code. You must also install the [python-matlab-bridge][3].

Make sure you can run the command `matlab` from a command window (or terminal).  If not, try the methods below:

__Unix-like-systems (including OSX):__

Add the following to your `~/.profile` or `~/.bash_profile`:

`export PATH="$PATH:<PATH_TO_MATLAB>"`

__Windows systems:__

From a CMD window:

`> setx PATH "%PATH%;C:\PATH_TO_MATLAB.EXE"`

Note, if you see the message:

"WARNING: The data being saved is truncated to 1024 characters"

It means your PATH variable is too long. You'll have to manually trim in in the Windows Environmental Variables editor.


Solutions of the exercises
------------

The solutions to the exercises [are available online](https://github.com/gpeyre/numerical-tours/tree/master/matlab/solutions).

[1]: https://github.com/gpeyre/numerical-tours/archive/master.zip
[2]: http://ipython.org/install.html
[3]: http://arokem.github.io/python-matlab-bridge/
