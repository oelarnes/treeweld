This folder contains software for calculating vertices and generating
images of conformally balanced trees of an arbitrary number of vertices.
This work is motivated by and developed as part of my PhD thesis.
Conformally balanced trees are, in short, a "true" representation of
a plane tree. Such trees are the pre-images of certain polynomials
characterized by the property of having all of its critical points
take on one of two values under the polynomial. Such polynomials are
topics of interest in algebraic geometry, but our own interest is
motivated by potential connections to the conformally invariant
processes of statistical physics. This arises from the connection
to the Brownian continuum random tree, and the problem of conformal
welding.

The program is divided into two scripts, laminations.py and weld.m. 

laminations.py

To generate a tree, run python laminations.py and follow the instructions
at the prompt. A reasonable set of values would be 20 edges, 1 tree, 10
subdivisions.

A lamination is encoded as a pairing of intervals. Consider
the rooted plane tree described by the following text diagram:

\_ _|_.
/   |

The '.' symbol represents the location of the root vertex. We can
encode the tree by proceeding counterclockwise around the outside
of the tree, assigning the numbers 1 to 14 to each side of each edge
in order. We then list the paired labels for each edge, generating
the following list for the tree above:

1 14
2 3
4 11
5 10
6 7
8 9
12 13

The script generates tables like the list above, with the
additional property that each pair is always a leaf of the 
sub-pairing consisting of the pair and all pairs below it.

Finally, given subdivisons, the program generates laminations
for the equivalent tree with vertices added in the middle of edges
corresponding to the specified number of subdivisions.


weld.m

Open Octave or MATLAB and run "weld" after generating the lamination tables
with laminations.py

weld.m performs the conformal mapping algorithm that generates
the tree figure in the plane. Points of the unit circle are generated
according to a beta distributions that ensures vertices will be
asymptotically evenly spaced in the image set. Then each leaf pair
(always the top pair of the lamination) is welded together via complex
function which has the following effect on the boundary:

1|
2|  |->  1|_2_
3|		 4| 3
4|

After all but the final pair are welded, the program maps back to the
unit circle and a final welding map is applied. The resulting line graph
obtained by connecting the images of successive points along the unit circle
is a tree in the complex plane, an approximation of the "true tree" with the
given combinatorics. For each possible combination, there exists precisely one
true tree with the specified normalization.

Finally, after obtaining a first approximation through the welding process, there
is an optional step to apply a vertex improvement scheme using Newton's method.
Currently this option is disabled, but using it we can obtain vertices accurate
to within machine accuracy for trees up to 2000 vertices.

tree.jpg

An example of a tree with many edges generated with weld.m

snowflake.jpg

A plot of the "b_1" values for each tree with 10 vertices. b_11 corresponds to the square
of the focus of the best-fit ellipse of a given tree. We have no theorem about the limiting distribution of b_1 for large trees, but this image and similar data suggest
uniform distribution near the origin, with either sharp or sigmoid decay at roughly 
|x| = .5

algorithms developed by:
Joel Barnes
Donald Marshall

coding:
Joel Barnes