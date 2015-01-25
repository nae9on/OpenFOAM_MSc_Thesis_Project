/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    
Description
    Calculates volume fractions from mass fractions
    It is possible to obtain the respective volume fractions of dry gas (without H2O)
    The "-dryGas" option requires a specie named "H2O"

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "hsReactionThermo.H"  //added
#include "RASModel.H" 
#include "edcChemistryModel.H"  
#include "thermoPhysicsTypes.H"
#include "reactingMixture.H" //added

using namespace Foam;



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    timeSelector::addOptions();    
    argList::validOptions.insert("dryGas", "");

#   include "setRootCase.H"

    //shall we remove water content? (=dryGas option set)
    bool dryGas = false;
    if (args.options().found("dryGas"))
    {
        dryGas = true;
   
    }
    
#   include "createTime.H"
#   include "createMesh.H"
    //access thermophysical model
    autoPtr<edcChemistryModel> edcChemistry
    (
        edcChemistryModel::New(mesh)
    );
    edcChemistryModel& chemistry = edcChemistry();
    hsReactionThermo& thermo = chemistry.thermo();


    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    //access specie thermo data 
    const PtrList<gasThermoPhysics> & specieThermo =
             dynamic_cast<const reactingMixture<gasThermoPhysics>&>  (thermo).speciesData();      
    
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );  
    

    
    //access U field (needed for turbulence model)
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    #   include "compressibleCreatePhi.H"
    
    //access turbulence model (needed for chemistry model)
    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo 
        )
    );



    // create new scalar field for mean molecular weight and initialize it to zero
    volScalarField Mmean
    (
        IOobject
        (
            "Mmean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Y[0]
    );
    
    //calculate Mmean from mass fractions
    //Mmean = (sum over all Yi/Mi )^(-1)
    Mmean=0.;

    forAll(Y,specieI)
    {
        forAll(Y[0].internalField(), cellI)
        {
            Mmean.internalField()[cellI]+=Y[specieI].internalField()[cellI]/specieThermo[specieI].W();
        }
        forAll(Y[0].boundaryField(), patchI)
        {
            
            forAll(Y[0].boundaryField()[patchI], faceI)
            {
                Mmean.boundaryField()[patchI][faceI]+=Y[specieI].boundaryField()[patchI][faceI]/specieThermo[specieI].W();
            }
        }        
            
    }


    forAll(Y[0].internalField(), cellI)
    {
        Mmean.internalField()[cellI]=1./Mmean.internalField()[cellI];
    }
    forAll(Y[0].boundaryField(), patchI)
    {
        
        forAll(Y[0].boundaryField()[patchI], faceI)
        {
            Mmean.boundaryField()[patchI][faceI]=1./Mmean.boundaryField()[patchI][faceI];
        }
    }

    //create new fields for volume fraction data
    PtrList<volScalarField> 	volY(Y.size());

    word prefix;
    if (dryGas)
    {
        prefix = "volDry";
    }
    else
    {
        prefix = "vol";        
    }

    forAll(Y, specieI)
    {

         volY.set
         (specieI,new volScalarField
             (
                 IOobject
                 (
                     prefix + Y[specieI].name(),
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 Y[0]*0.
              )
          );
     }    

    
    //prepare for dryGas option

    // find H2O index
    label indexH2O = -1;
    forAll(Y,specieI)
    {
        if (Y[specieI].name()== "H2O" )
        {
            indexH2O = specieI;
            break;
        }
    }
    if(dryGas && (indexH2O==-1))
    {
             FatalErrorIn("massToVolumeFraction:")
                 << "Could not find specie with name 'H2O' which is required for option '-dryGas' "
                 << exit(FatalError);
    }
    volScalarField volH2O = Mmean/specieThermo[indexH2O].W() * Y[indexH2O];
    //prevent division by zero
    volH2O=min(volH2O, 0.999999999); 
    
    //calculate and write volume fraction data
    //correct for dryGas, if neccessary
    forAll(Y,specieI)
    {
        
        forAll(Y[0].internalField(), cellI)
        {
            volY[specieI].internalField()[cellI]=Mmean.internalField()[cellI]/specieThermo[specieI].W() * Y[specieI].internalField()[cellI];
            if(dryGas)
            {
                volY[specieI].internalField()[cellI] /= (1-volH2O.internalField()[cellI]);
            }            
        }
        forAll(Y[0].boundaryField(), patchI)
        {
            
            forAll(Y[0].boundaryField()[patchI], faceI)
            {
                volY[specieI].boundaryField()[patchI][faceI]=Mmean.boundaryField()[patchI][faceI]/specieThermo[specieI].W() * Y[specieI].boundaryField()[patchI][faceI];
            if(dryGas)
            {
                volY[specieI].boundaryField()[patchI][faceI] /= (1-volH2O.boundaryField()[patchI][faceI]);
            }              

            }
        }
        
        // write volScalarFields
        if (!dryGas)
        {
           volY[specieI].write();        
        }
        else // do not write volDryH2O, since it makes no sense
        {
            if (specieI!=indexH2O)
               volY[specieI].write();             
        }
            
    }
    
    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //

