---
layout: page
title: "Julia's tours"
description: "Numerical Tours in Julia"
header-img: "img/hokusai-8.jpg"
---

These are the [Julia](http://julialang.org/) tours, that can be browsed as HTML pages, but can also be downloaded as [Jupyter notebooks](http://jupyter.org/). Please read the [installation page](../installation_julia/) for more information about how to run these tours.


{% for categ in site.data.tours_julia %}

{{ categ.name }}      {#{{ categ.short }}}
----------------

<ul>
{% for tour in categ.list %}
	{% if tour.format == "html" %}
		<li> <a href="{{ tour.rep }}"> {{ tour.name }} </a> </li>
	{% else %}
		<li> <a href="http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/julia/{{ tour.rep }}.ipynb"> {{ tour.name }} </a> </li>
	{% endif %}
{% endfor %}
</ul>

{% endfor %}
