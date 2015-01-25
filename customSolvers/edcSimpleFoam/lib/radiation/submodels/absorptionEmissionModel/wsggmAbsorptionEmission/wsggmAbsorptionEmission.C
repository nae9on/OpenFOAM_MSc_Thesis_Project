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

\*---------------------------------------------------------------------------*/

#include "wsggmAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(wsggmAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            wsggmAbsorptionEmission,
            dictionary
        );
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmission::wsggmAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("e", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    emissivityCoeffs_(coeffsDict_.lookup("emissivityCoeffs")),
    fittingFactors_(coeffsDict_.lookup("fittingFactors")),
    pathLength_(coeffsDict_.lookup("pathLength"))    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggmAbsorptionEmission::~wsggmAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::aCont(const label bandI) const
{

    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "a",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );    

       
    const hsReactionThermo& thermo= mesh_.lookupObject<hsReactionThermo>("thermophysicalProperties");

    //access species mass fractions
    const PtrList<volScalarField>& Y = thermo.composition().Y();
    //access specie thermo data 
    const PtrList<gasThermoPhysics> & specieThermo =
        dynamic_cast<const reactingMixture<gasThermoPhysics>&>  (thermo).speciesData();
    // get index of CO2 in mixture
    label indexCO2= dynamic_cast<const reactingMixture<gasThermoPhysics>&> (thermo).species()["CO2"];    
    // get index of H2O in mixture
    label indexH2O= dynamic_cast<const reactingMixture<gasThermoPhysics>&> (thermo).species()["H2O"];    


    volScalarField emissivity(        
            IOobject
            (
                "emissivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    volScalarField pressurePathLength(        
            IOobject
            (
                "pressurePathLength",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimPressure*dimLength,0.0)
        );
    volScalarField weightingFactor(        
            IOobject
            (
                "weightingFactor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    const scalar wCO2=44;
    const scalar wH2O=18;

    volScalarField wMean(        
            IOobject
            (
                "wMean",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0),0.0) // kg/kmol
        );

    forAll(Y,specieI)
    {
        wMean+=Y[specieI]/specieThermo[specieI].W();
    }

    wMean=1/wMean;

    pressurePathLength=wMean*thermo.p()*(Y[indexCO2]/wCO2+Y[indexH2O]/wH2O)*pathLength_;

    volScalarField limitedDimlessTemperature = min(thermo.T()/dimensionedScalar("unity",dimTemperature, 1.0), 2400.0);

    forAll(emissivityCoeffs_,i)
    {
        weightingFactor = 0.0;
        for(label j=1; j<fittingFactors_.size(); j++)
        {
            weightingFactor += fittingFactors_[i][j]*pow(limitedDimlessTemperature, (j-1));
        }
        emissivity += weightingFactor * (1 - exp( (-1) *emissivityCoeffs_[i] * pathLength_.value()));
    }
//    a_= (-1)* log(1-emissivity) / pathLength_;
    emissivity=min(emissivity,0.9999);
    ta()= (-1)* log(1-emissivity) / pathLength_;

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::eCont(const label bandI) const
{
  
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "e",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    const hsReactionThermo& thermo= mesh_.lookupObject<hsReactionThermo>("thermophysicalProperties");

    //access species mass fractions
    const PtrList<volScalarField>& Y = thermo.composition().Y();
    //access specie thermo data 
    const PtrList<gasThermoPhysics> & specieThermo =
        dynamic_cast<const reactingMixture<gasThermoPhysics>&>  (thermo).speciesData();
    // get index of CO2 in mixture
    label indexCO2= dynamic_cast<const reactingMixture<gasThermoPhysics>&> (thermo).species()["CO2"];    
    // get index of H2O in mixture
    label indexH2O= dynamic_cast<const reactingMixture<gasThermoPhysics>&> (thermo).species()["H2O"];    


    volScalarField emissivity(        
            IOobject
            (
                "emissivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    volScalarField pressurePathLength(        
            IOobject
            (
                "pressurePathLength",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimPressure*dimLength,0.0)
        );
    volScalarField weightingFactor(        
            IOobject
            (
                "weightingFactor",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimless,0.0)
        );
    const scalar wCO2=44;
    const scalar wH2O=18;

    volScalarField wMean(        
            IOobject
            (
                "wMean",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0),0.0) // kg/kmol
        );

    forAll(Y,specieI)
    {
        wMean+=Y[specieI]/specieThermo[specieI].W();
    }

    wMean=1/wMean;

    pressurePathLength=wMean*thermo.p()*(Y[indexCO2]/wCO2+Y[indexH2O]/wH2O)*pathLength_;

    volScalarField limitedDimlessTemperature = min(thermo.T()/dimensionedScalar("unity",dimTemperature, 1.0), 2400.0);

    forAll(emissivityCoeffs_,i)
    {
        weightingFactor = 0.0;
        for(label j=1; j<fittingFactors_.size(); j++)
        {
            weightingFactor += fittingFactors_[i][j]*pow(limitedDimlessTemperature, (j-1));
        }
        emissivity += weightingFactor * (1 - exp( (-1) *emissivityCoeffs_[i] * pathLength_.value()));
    }
//    a_= (-1)* log(1-emissivity) / pathLength_;
    emissivity=min(emissivity,0.9999);
    te()= (-1)* log(1-emissivity) / pathLength_;

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggmAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE=E_;
 

    return tE;

}


// ************************************************************************* //
