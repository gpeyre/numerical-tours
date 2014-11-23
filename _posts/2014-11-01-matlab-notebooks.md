---
layout:     post
title:      "Matlab notebooks"
subtitle:   "using iPython notebooks"
date:       2014-11-01 12:00:00
author:     "Gabriel Peyré"
header-img: "/img/hokusai-12.jpg"
---

Thanks to [Steven Silvester](https://github.com/blink1073), most of the [Matlab tours]({{ site.baseurl }}/matlab/) are now available as iPython notebooks. This means that:

* The corresponding tours are now directly rendered online using [nbviewer](http://nbviewer.ipython.org/), which means that HTML conversion is not anymore needed.
* The user can now download the corresponding .ipynb file and run it locally, possibly modifying and completing its content.

This requires that you install Python and iPython (I recommend for instance the [Anaconda](http://continuum.io/downloads) distribution). You also need to install the [pymatbridge](https://pypi.python.org/pypi/pymatbridge) Python module in order to run Matlab code from the notebook.

The command 

> %load_ext pymatbridge

will start Matlab in the background, and then you simply need to put the keyword 

> %%matlab

if you want to add a new cell to the notebook with your own Matlab code.   