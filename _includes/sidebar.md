{% capture lvl %}{{ page.url | append:'index.html' | split:'/' | size }}{% endcapture %}
{% capture relative %}{% for i in (3..lvl) %}../{% endfor %}{% endcapture %}

[home]({{ relative }}index.html)

**[install]({{ relative }}doc/install.html)**

[documentation]({{ relative }}doc/index.html)

[papers]({{ relative }}papers/papers.html)
