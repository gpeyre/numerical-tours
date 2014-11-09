---
layout: page
title: "Matlab's tours"
description: "Numerical Tours in Matlab"
header-img: "/img/hokusai-6.jpg"
---

<p align="center">
<br/>
<br/>
{% for categ in site.data.tours %}
<a href="#{{ categ.short }}"> {{ categ.name }} </a> <br/>
{% endfor %}
</p>

<br/><br/>


{% for categ in site.data.tours %}

{{ categ.name }}      {#{{ categ.short }}}
----------------

<ul>
{% for tour in categ.list %}
	<li> <a href="http://nbviewer.ipython.org/github/gpeyre/numerical-tours/blob/master/matlab/{{ tour.rep }}.ipynb"> {{ tour.name }} </a> </li>
{% endfor %}
</ul>

{% endfor %}
