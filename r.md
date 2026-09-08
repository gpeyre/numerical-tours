---
layout: page
title: "R archive"
description: "Numerical Tours in R"
header-img: "img/hokusai-11.jpg"
---

<div class="archive-notice">This is the R collection. For complete Python walkthroughs, visit the <a href="{{ '/python/' | relative_url }}">Python tours</a>. <a href="{{ '/archive/' | relative_url }}">All language archives →</a></div>


These are the [R](https://www.r-project.org/) tours, that can be browsed as HTML pages, but can also be downloaded as [Jupyter notebooks](https://jupyter.org/). Please read the [installation page](../installation_r/) for more information about how to run these tours.


{% for categ in site.data.tours_r %}

{{ categ.name }}      {#{{ categ.short }}}
----------------

<ul>
{% for tour in categ.list %}
	{% if tour.format == "html" %}
		<li> <a href="{{ tour.rep }}"> {{ tour.name }} </a> </li>
	{% else %}
		<li> <a href="https://nbviewer.org/github/gpeyre/numerical-tours/blob/master/r/{{ tour.rep }}.ipynb"> {{ tour.name }} </a> </li>
	{% endif %}
{% endfor %}
</ul>

{% endfor %}
