Draw polyhedra using OpenGL.

OpenGL 2.2 ("modern" version forthcoming).

Uses C++11 features.

Usage eg. with 3 args:

./poly --zrot 4 2 3 3 2 3 5

Displays polyhedron based on 4 2 3 Schwarz triangle, with additional compounds according to 2 3 5 symmetry. With "rotations only" ('x' key) and compound colours ('f' key), shows compound of 5 tetrahedra.

Various keyboard commands:

]: move to next position in Schwarz triangle
[: move to previous position in Schwarz triangle
space: toggle continuous movement around Schwarz triangle
c: toggle compounding
z: toggle z rotation for compounds
f: toggle face/compound colour mode
d: change dual mode
k: change dual style
r: toggle line drawing mode
<: move in
>: move out
<arrow keys>: rotate
p: stop rotation
x: toggle rotations only for compounding


