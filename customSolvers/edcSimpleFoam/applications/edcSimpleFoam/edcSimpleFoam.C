/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    rhoSimpleFoam

Description
    Steady-state SIMPLE solver for reacting laminar or turbulent RANS flow of
    compressible fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "basicPsiThermo.H"   //removed
#include "RASModel.H"
//#include "hCombustionThermo.H"  //added
#include "hsReactionThermo.H"  //added
#include "edcChemistryModel.H"  //added
#include "multivariateScheme.H" //added
#include "thermoPhysicsTypes.H" //added
#include "Reaction.H" //added
#include "reactingMixture.H" //added
#include "radiationModel.H" //added
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createRadiationModel.H"     //added
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        

        #include "readSIMPLEControls.H"
        #include "initConvergenceCheck.H"

        p.storePrevIter();
        rho.storePrevIter();

        #include "chemistry.H" //added

        // Pressure-velocity SIMPLE corrector
        {
            
            #include "UEqn.H"
            #include "pEqn.H"

            #include "multivariateConvection.H"            
            #include "hsEqn.H"
            #include "YEqn.H"   //added        

        }



        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        #include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
