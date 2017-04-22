This is a post processing application that calculates the terms in the X-component of the integral form of the vorticity transport equation. This application saves the computed fields and writes it in the respective time directories. The following quantities are written:-
```
1) vorticity and magVorticity
2) vortConvection -- convection of vorticity term
3) vortDiffusion  -- Diffusion of vorticity term
4) shearFlux      -- the shear layer flux of vorticity
5) yTilting       -- term signifying the y tilt of vorticity
6) zTilting       -- term signifying the z tilt of vorticity 
```
An integral analysis of the vorticity transport equation can then be performed by
sampling these fields over definitive contours and taking the surface area integral.

# Installation #

Source the openfoam environment (bashrc or cshrc file) to load the proper environment variables. Place this folder in the $WM_PROJECT_USER_DIR/utilities/ directory and run wmake from the source directory. The application will be installed in the directory set up by the environment variable $FOAM_USER_APPBIN

# Usage #
For usage help, run vorticityTransportTerms -help