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
    edcPisoFoam

Description
    -Transient solver for gas combustion using edcChemistryModel and thus providing 
     access to several flavors of the Eddy Dissipation Concept (EDC) Model by B.Magnussen
    
    -Time step adaption is supported


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "turbulenceModel.H"
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
    #include "readGravitationalAcceleration.H"    
    #include "initContinuityErrs.H"

    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        #include "readTimeControls.H"
        #include "readPISOControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"       

        #include "initConvergenceCheck.H"
        
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;        

        #include "chemistry.H" //added
        rho = thermo.rho();
        #include "rhoEqn.H"

        for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
        {
            #include "UEqn.H"
            #include "multivariateConvection.H" //has to be called before YEqn.H and hsEqn.H 
            #include "YEqn.H"
            #include "hsEqn.H"

            // --- PISO loop
            for (int corr=1; corr<=nCorr; corr++)
            {
                #include "pEqn.H"
            }
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
