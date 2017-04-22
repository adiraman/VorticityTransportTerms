//Define velocity components:-
volScalarField  UX=U.component(0);
volScalarField  UY=U.component(1);
volScalarField  UZ=U.component(2);

//Define velocity gradient components:-
volScalarField  gradUXX=gradU.component(0);
volScalarField  gradUYY=gradU.component(4);
volScalarField  gradUZZ=gradU.component(8);

volScalarField  gradUXY=gradU.component(1);
volScalarField  gradUYZ=gradU.component(5);
volScalarField  gradUZX=gradU.component(6);

volScalarField  gradUYX=gradU.component(3);
volScalarField  gradUZY=gradU.component(7);
volScalarField  gradUXZ=gradU.component(2);

//Define vorticity components:-
volScalarField  vortX=U.component(0);
volScalarField  vortY=U.component(1);
volScalarField  vortZ=U.component(2);

//Define vorticity gradient components:-
volScalarField  gradVortXX=gradVort.component(0);
volScalarField  gradVortYY=gradVort.component(4);
volScalarField  gradVortZZ=gradVort.component(8);

volScalarField  gradVortXY=gradVort.component(1);
volScalarField  gradVortYZ=gradVort.component(5);
volScalarField  gradVortZX=gradVort.component(6);

volScalarField  gradVortYX=gradVort.component(3);
volScalarField  gradVortZY=gradVort.component(7);
volScalarField  gradVortXZ=gradVort.component(2);

