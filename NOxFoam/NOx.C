/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvIOoptionList.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "CourantNo.H"

    #include "createSpecies.H"

    #include "preProcessData.H"

    if (flagPrompt.value()) 
    {
    #include "createSourcePromptNO.H"
    }

    #include "updateSourceThermalNO.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(rhoAvg,NO)
              + fvm::div(fvc::interpolate(rhoAvg)*phi, NO)
              - fvm::laplacian(rhoDTeff, NO)
             ==
                flagThermal*sourceThermalNO + flagPrompt*sourcePromptNO
                //fvOptions(NO) Source(NO) sourceThermalNO sourcePromptNO 
            );
        }
        
        // updating total NO ppm and concentration
        ppmNO = NO*pow(10,6);
        concNO = NO*p/(Rgas*T);
        if(runTime.outputTime()) concNO.write();

        Info<< "Source for ThermalNO updated" << endl;
        #include "updateSourceThermalNO.H"
        
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
