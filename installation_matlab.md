---
layout: page
title: "Installation"
description: "for Matlab"
header-img: "/img/hokusai-5.jpg"
---

Installation
------------

You need to download [numerical_tours][1]
and install the IPython [notebook][2] to run the code.

You must also install the [python-matlab-bridge][3].

Make sure you can run the command `matlab` from a command window
(or terminal).  If not, try the methods below:

__Unix-like-systems (including OSX):__

Add the following to your `~/.profile`: 

`export PATH="$PATH:<PATH_TO_MATLAB>"`

__Windows systems:__

From a CMD window: 

`> setx PATH "%PATH%;C:\PATH_TO_MATLAB.EXE"`

Note, if you see the message: 

"WARNING: The data being saved is truncated to 1024 characters" 

It means your PATH variable is too long. You'll have to manually trim in in the Windows Environmental Variables editor.

[1]: https://github.com/gpeyre/numerical-tours/archive/master.zip
[2]: http://ipython.org/install.html
[3]: http://arokem.github.io/python-matlab-bridge/