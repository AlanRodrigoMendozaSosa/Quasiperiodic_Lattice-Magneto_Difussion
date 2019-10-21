# Quasiperiodic Lattice

The main objective of the code is to study the diffusion of an electron in a two-dimensional plate doped with quantum antidots in the presence of an external magnetic field. The centers of the quantum antidots correspond to the positions of a quasiperiodic lattice.

To simulate that, we use a generalization of the "Cut-and-Project" method to obtain the positions of the 2D quasiperiodic lattice.

The central idea is to generate a cylinder with an irrational inclination with respect to the XY plane inside a unit cube centered on the origin with periodic border conditions. That cylinder will simulate our quantum antidots when we work only in a 2D projection of the 3D space.

Now we put a point inside the cube, that point can only move in planes orthogonal to the cylinder, following circular paths. That point will simulate our electron moving in a 2D plate in the presence of an external magnetic field.

If the particle collides with the cylinder, then the particle will change his direction following a specular-like reflexion.

An extended explanation about the algorithm used can be found in the author thesis named ["Movimiento circular en ambientes cuasi-periódicos: magnetorresistencia en cuasicristales"](http://oreon.dgbiblio.unam.mx/F/JA719DL9AX5IAJIN99N6FJLDSVNI8IQU6FYCSM4VV12VL8C94X-41442?func=find-b&request=Alan+Rodrigo+Mendoza+Sosa&find_code=WAT&adjacent=N&local_base=TES01&x=46&y=11&filter_code_2=WYR&filter_request_2=&filter_code_3=WYR&filter_request_3=)
