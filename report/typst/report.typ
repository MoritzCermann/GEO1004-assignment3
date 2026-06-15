#set page(
  numbering: "1"
)

#set page(columns: 2)

#show figure.caption: set text(size: 10pt, style: "italic")

#set text(
  size: 10pt,
  lang: "en", 
  region: "gb"
)
#set par(justify: true)

#show title: set text(size: 17pt)
#show title: set block(below: 1.2em)

#set heading(numbering: "1.")
#show heading.where(level: 1): set text(size: 14pt)

#place(top + center,
  float: true, 
  scope: "parent",
  clearance: 3em,)[

#align(center)[
  #title[GEO1004 Assignment 3: Simple Polyhedron Processing]
  _#datetime.today().display("[day padding:none] [month repr:long] [year]")_
  #grid(columns: (1fr, 1fr, 1fr),
  [Ruben Vons \ 6285880],
  [Artemi Kurski \ 6195784],
  [Moritz Cermann \ 6562906]) 
  ]
]



#include "content.typ"