## Draw polyhedra using OpenGL.

Now superceded by WebGL/THREE.js version.

Uses old-style fixed pipeline OpenGL. 

Uses C++11 features.

Usage eg. with 3 args:

./poly 2 3 5

Displays polyhedron based on 2 3 5 Schwarz triangle (ie. basic icosahedral symmetry).

./poly --zrot 4 2 3 3 2 3 5

Basic polyhedron 2 3 3/2, with compounding according to 2 3 5 symmetry.

Various keyboard commands:

* \]: move to next position in Schwarz triangle
* \[: move to previous position in Schwarz triangle
* space: toggle continuous movement around Schwarz triangle
* c: toggle compounding
* z: toggle z rotation for compounds
* f: toggle face/compound colour mode
* d: change dual mode
* k: change dual style
* r: toggle line drawing mode
* <: move in
* \>: move out
* \[arrow keys\]: rotate
* p: stop rotation
* x: toggle rotations only for compounding


