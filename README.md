# Globetrotter

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11668651.svg)](https://doi.org/10.5281/zenodo.11668651)

Globetrotter calculates weighed cumulative distance of each grid cell to the 
nearest point in a set of input points or grid cells based on a weight grid
(i.e. unit travel cost grid). 

Multiple input locations are processed at once by using a two-phase concurrent
implementation of [Dijkstra's algorithm][1]. In the first phase, the algorithm
divides the grid into multiple mutually exclusive blocks which are processed
individually in parallel. In the second phase, block results are merged and
corrected by re-considering cells located at the boundaries of the blocks.

In addition to the minimum weighted distance grid, the algorithm also
returns a grid of nearest reference location ids, which enables the catchments
of the reference locations to be identified easily.

Both Euclidean and spherical coordinate systems are supported. It is also 
possible to roll in x and y axes, i.e. at the International Date Line and poles.

The algorithm is designed to be fast and memory efficient. Therefore, it
is suitable not only for moderate sized areas, but also global studies
(e.g. global accessibility map).

[1]: https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
