/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    vorticity

Group
    grpPostProcessingUtilities

Description
    Calculates and writes the vorticity of velocity field U.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject vorticityHeader
    (
        "vorticity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (Uheader.typeHeaderOk<volVectorField>(true))
    {
        Info<< "    Reading U" << endl;
        Info<< " " << endl;
        volVectorField U(Uheader, mesh);

        volVectorField* vortPointer;

        if (vorticityHeader.typeHeaderOk<volVectorField>(true))
        {
            Info<< "    Reading vorticity field" << endl;
            Info<< " " << endl;
            vortPointer = new volVectorField
            (
                IOobject
                (
                    "vorticity",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            );
        }
        else
        {
            Info<< "    Calculating vorticity field" << endl;
            vortPointer = new volVectorField
            (
                IOobject
                (
                    "vorticity",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::curl(U)
            );

            volScalarField magVorticity
            (
                IOobject
                (
                    "magVorticity",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                mag(*vortPointer)
            );
            Info<< "vorticity max/min : "
                << max(magVorticity).value() << " "
                << min(magVorticity).value() << endl;
            Info<< " " << endl;

            if (writeResults)
            {
                magVorticity.write();
            }
        }

        Info<< "    Calculating vorticity transport terms..." << endl;
        Info<< " " << endl;
        volTensorField gradU = fvc::grad(U);
        volTensorField gradVort = fvc::grad(*vortPointer);
        #include "gradients.H"

        Info << "    Calculating Vortex Convection Term" << endl;
        volScalarField vortConvection
        (
            IOobject
            (
                "vortConvection",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            gradVortXX*UX
        );

        Info<< "vortex convection max/min: "
            << max(vortConvection).value() << " "
            << min(vortConvection).value() << endl;
        Info<<" " << endl;

        Info << "    Calculating Vorticity source from Y-Tilting" << endl;
        volScalarField yTilting
        (
            IOobject
            (
                "yTilting",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            vortY*gradUXY
        );

        Info<< "vortex Y tilting max/min: "
            << max(yTilting).value() << " "
            << min(yTilting).value() << endl;
        Info<< " " << endl;

        Info<< "    Calculating Vorticity source from Z-Tilting" << endl;
        volScalarField zTilting
        (
            IOobject
            (
                "zTilting",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            vortZ*gradUXZ
        );

        Info<< "vortex Z tilting max/min: "
            << max(zTilting).value() << " "
            << min(zTilting).value() << endl;
        Info<<" " << endl;

        // This is a round about fix to solve the dimension mismatch error!
        dimensionedScalar fix("fix", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1.0);

        Info<< "    Calculating the Shear Layer flux to vorticity" << endl;
        volScalarField shearFlux
        (
            IOobject
            (
                "shearFlux",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            UY*gradVortXY + UZ*gradVortXZ - (vortX*gradUXX/fix)
        );

        Info<< "vortex shear Flux max/min: "
            << max(shearFlux).value() << " "
            << min(shearFlux).value() << endl;
        Info<<" " << endl;

        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dimensionedScalar nu
        (
            transportProperties.lookup("nu")
        );

        Info<< "    Calculating the Diffusion of Vorticity" << endl;
        volScalarField vortDiffusion
        (
            IOobject
            (
                "vortDiffusion",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            nu*fvc::laplacian(vortX)
        );

        Info<< "vortex diffusion max/min: "
            << max(vortDiffusion).value() << " "
            << min(vortDiffusion).value() << endl;

        if (writeResults)
        {
            vortPointer->write();
            vortConvection.write();
            yTilting.write();
            zTilting.write();
            shearFlux.write();
            vortDiffusion.write();
        }

        delete vortPointer;
    }
    else
    {
        Info<< "    No U" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
