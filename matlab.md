---
layout: page
title: "Matlab's tours"
description: "Numerical Tours in Matlab"
header-img: "img/hokusai-6.jpg"
---

These are the [Matlab](http://www.mathworks.fr/products/matlab/) tours, that can be browsed as HTML pages, but can also be downloaded as iPython notebooks. Please read the [installation page](../installation_matlab/) for more information about how to run these tours. A lot of Matlab tours are also compatible with [Scilab](http://www.scilab.org/) and with [GNU Octave](https://www.gnu.org/software/octave/).

<p align="center">
<br/>
<br/>
{% for categ in site.data.tours_matlab %}
<a href="#{{ categ.short }}"> {{ categ.name }} </a> <br/>
{% endfor %}
</p>

<br/><br/>


{% for categ in site.data.tours_matlab %}

{{ categ.name }}      {#{{ categ.short }}}
----------------

<ul>
{% for tour in categ.list %}
	<li>
	{% if tour.format == "html" %}
		<a href="{{ tour.rep }}"> {{ tour.name }} </a>
	{% else %}
		<a href="http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/matlab/{{ tour.rep }}.ipynb"> {{ tour.name }} </a>
	{% endif %}
	&nbsp;&nbsp;(<a href="http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/matlab/{{ tour.rep }}.ipynb">ipynb</a>|<a href="{{ tour.rep }}">web</a>)
	</li>
{% endfor %}
</ul>

{% endfor %}
