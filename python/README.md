About 
======

This directory stores the iPython notebook (.ipynb) of the "Numerical Tours" (www.numerical-tours.com). 

Usage
======

By clicking on a .ipynb file, you will access a static, compiled, version of the tour. You can download the file in order to run dynamically the notebook (using "ipython notebook --pylab inline" to start ipython notebooks), and modify its content. For this to work, you will also need to download the nt_toolbox.py file. 

How to contribute
======

This is a work in progress to port all the Numerical Tours to Python. All the un-ported tours are in the directory "todo/". These consist in raw .ipynb files that have been exported from the corresponding .m file using the tool nb_converter.py provided by Steven Silvester.

In order to help in this porting task, you need to select your favorite un-ported tour, edit it so that the notebook works correctly, and then move the corrected .ipynb file from python/todo/ to python/.

This might require porting some functions from the matlab toolboxes (toolbox_signal/, toolbox_general/, toolbox_graph/) to python. The newly created functions should be added to the file nt_toolbox.py.

======

Copyright (c) 2014 Gabriel Peyre