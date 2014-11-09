---
layout: page
title: "Matlab's tours"
description: "Numerical Tours in Matlab"
header-img: "/img/hokusai-6.jpg"
---

These are the [Matlab](http://www.mathworks.fr/products/matlab/) tours, that can be browsed as HTML pages, but can also be downloaded as iPython notebooks. Please read the [installation page](../installation_matlab/) for more information about how to run these tours.

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
	{% if tour.format == "html" %}
		<li> <a href="{{ tour.rep }}"> {{ tour.name }} </a> </li>
	{% else %}
		<li> <a href="http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/matlab/{{ tour.rep }}.ipynb"> {{ tour.name }} </a> </li>
	{% endif %}
{% endfor %}
</ul>

{% endfor %}
