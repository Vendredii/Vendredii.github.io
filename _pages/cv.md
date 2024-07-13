---
layout: archive
title: "CV"
permalink: /cv/
author_profile: true
redirect_from:
  - /resume
---

{% include base_path %}

Education
======
* Ph.D student in Ecology, Utrecht University, 2026 (expected)
* M.S. in Ecology, Nanjing Forestry University, 2022
* B.S. in Forestry, Southwest Forestry University, 2019
  
Skills
======
* Programming
  * R
  * Python
* GIS analysis
  * ArcGIS Pro
  * ESRI
  * Google Earth Engine

Publications
======
  <ul>{% for post in site.publications reversed %}
    {% include archive-single-cv.html %}
  {% endfor %}</ul>
  
Talks
======
  <ul>{% for post in site.talks reversed %}
    {% include archive-single-talk-cv.html  %}
  {% endfor %}</ul>
  
Teaching
======
  <ul>{% for post in site.teaching reversed %}
    {% include archive-single-cv.html %}
  {% endfor %}</ul>
  
Service and leadership
======
* Organizing cannoing with my colleagues?
