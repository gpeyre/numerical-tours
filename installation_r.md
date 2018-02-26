---
layout: page
title: "Installation"
description: "for R"
header-img: "img/hokusai-5.jpg"
---

Installation of IRkernel
------------

The tours in R are Jupyter notebooks. To run them (which we recommend), you need to install IRkernel through https://github.com/IRkernel/IRkernel.


Installation of the libraries
------------

The tours make use of different external libraries such as:

* imager
* pracma
* Matrix

Each time you encounter a new library in a numerical tour, you will need to install it.\
In order to install each of them on the R kernel for Jupyter, you need to run the following line on an R terminal:\
```
install.packages("name_of_library","directory_of_IRkernel_libraries", dep=TRUE)
```

The directory_of_IRkernel_libraries is the place of the subfolder R\library inside the installation folder of Anaconda3 on your computer.
Example of the installation of imager:
```
install.packages("imager","C:\Program Files\Anaconda3\R\library", dep=TRUE)
```

You can use from the R prompt the command home to locate this directory:
```
R.home()
```

Installation of the toolboxes packages
------------

Each tour makes use of a toolbox package that needs to be downloaded as a .zip file, and then unzipped within your working directory (from which you will run the notebooks): (https://github.com/gpeyre/numerical-tours/tree/master/r/nt_toolbox).

Alternatively, you can download the whole numerical_tours: (https://github.com/gpeyre/numerical-tours/tree/master/r/).


Solutions of the exercises
------------

The solutions to the exercises [are available online](https://github.com/gpeyre/numerical-tours/tree/master/r/nt_solutions).
