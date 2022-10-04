---
layout: page
title: "Python's tours"
description: "Numerical Tours in Python"
header-img: "img/hokusai-2.jpg"
---

These are the [Python](https://www.python.org/) tours, that can be browsed as HTML pages, but can also be downloaded as [Jupyter notebooks](http://jupyter.org/). Please read the [installation page](../installation_python/) for more information about how to run these tours.


{% for categ in site.data.tours_python %}

{{ categ.name }}      {#{{ categ.short }}}
----------------

<ul>
{% for tour in categ.list %}
	{% if tour.format == "html" %}
		<li> <a href="{{ tour.rep }}"> {{ tour.name }} </a> </li>
	{% else %}
		<li> <a href="http://nbviewer.org/github/gpeyre/numerical-tours/blob/master/python/{{ tour.rep }}.ipynb"> {{ tour.name }} </a> </li>
	{% endif %}
{% endfor %}
</ul>

{% endfor %}
