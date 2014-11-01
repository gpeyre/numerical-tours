---
layout: page
title: "Python's tours"
description: "Numerical Tours in Python"
header-img: "/img/hokusai-2.jpg"
---

These are the Python tours, that can be browsed as HTML pages, but can also be downloaded as iPython notebooks. Please read the [installation page](../installation_python/) for more information about how to run these tours.

Note that it is a work in progress to port all the Numerical Tours to Python. Help is wellcome, please refer to the GitHub repository for how to proceed.


{% for categ in site.data.tours_matlab %}

{{ categ.name }}      {#{{ categ.short }}}
----------------

<ul>
{% for tour in categ.list %}
	{% if tour.format == "html" %}
		<li> <a href="{{ tour.rep }}"> {{ tour.name }} </a> </li>
	{% else %}
		<li> <a href="http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/python/{{ tour.rep }}.ipynb"> {{ tour.name }} </a> </li>
	{% endif %}
{% endfor %}
</ul>

{% endfor %}




