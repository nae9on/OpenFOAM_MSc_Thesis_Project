/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "edm.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edm<CompType, ThermoType>::edm
(
    const fvMesh& mesh,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    CompType(mesh, thermoTypeName),

    Y_(this->thermo().composition().Y()),

    reactions_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>(this->thermo())
    ),
    specieThermo_
    (
        dynamic_cast<const reactingMixture<ThermoType>&>
            (this->thermo()).speciesData()
    ),

    nSpecie_(Y_.size()),
    nReaction_(reactions_.size()),


    RR_(nSpecie_)

{
    // create the fields for the chemistry sources
    forAll(RR_, fieldI)
    {
        RR_.set
        (
            fieldI,
            new scalarField(mesh.nCells(), 0.0)
        );
    }

    Info<< "edm: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edm<CompType, ThermoType>::~edm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edm<CompType, ThermoType>::tc() const
{
    FatalError << "edm::tc called, but not implemented!!"  << exit(FatalError);

    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edm<CompType, ThermoType>::dQ() const
{
    FatalError << "edm::dQ called, but not implemented!!"  << exit(FatalError);

    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edm<CompType, ThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        scalarField& Sh = tSh();
        forAll(Y_, i)
        {
            forAll(Sh, cellI)
            {
                scalar hi = specieThermo_[i].Hc() ;
                Sh[cellI] -= hi*RR_[i][cellI];
            }
        }
    }


    return tSh;
}

template<class CompType, class ThermoType>
Foam::scalar Foam::edm<CompType, ThermoType>::solve
(
    const scalar t0,
    const scalar deltaT
)
{   
    if(!this->chemistry_)
    {
        Info << "Chemistry disabled in const/chemistryProperties!!" << endl;
        return 0;
    }
    Info << "Selected chemistry model is Edm/Kinetic" << endl;
    //calculate EDM rates
    
    //model constants
    const scalar AEdm = 4.0;
    const scalar BEdm = 0.5;
    
    volScalarField& rho = this->thermo().rho();
    const volScalarField& k = this->db().objectRegistry::lookupObject<volScalarField>("k");         
    const volScalarField& epsilon = this->db().objectRegistry::lookupObject<volScalarField>("epsilon");             
    
    scalarField YMin(rho.size(),1.);
    scalarField edmFactor(rho.size(),0.);

  

    forAll(reactions_ , reactionI)
    {
        // minimum mass fraction must be smaller or equalthan "1." Set YMin to this max. value first
        YMin=1.;
        // determine minimum mass fraction from reactants
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        {
            label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;
            scalar stoichCoeff = reactions_[reactionI].lhs()[lhsSpecieI].stoichCoeff;
            if (!(stoichCoeff == 0))
            {
                YMin =  min(YMin, ( Y_[specieI].internalField() / stoichCoeff / specieThermo_[specieI].W() )  );
            }
        }

        // determine minimum mass fraction from products
        scalarField YPMin(rho.size(),0.);
        forAll(reactions_[reactionI].rhs() , rhsSpecieI)
        {
            label specieI = reactions_[reactionI].rhs()[rhsSpecieI].index;
            scalar stoichCoeff = reactions_[reactionI].rhs()[rhsSpecieI].stoichCoeff;
            if (!(stoichCoeff == 0))
            {
                YPMin+=  Y_[specieI].internalField() / stoichCoeff / specieThermo_[specieI].W() ;
            }
        }
        YMin =  min(YMin, YPMin*BEdm );
        
        //calculate reaction rate according to EDM model
        edmFactor=  AEdm * rho.internalField() * epsilon.internalField()/k.internalField() * YMin;

        //create a source term according to EDM reaction rate for each species
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        {
                label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;
                scalar stoichCoeff = reactions_[reactionI].lhs()[lhsSpecieI].stoichCoeff;
                RR_[specieI]-= stoichCoeff  * specieThermo_[specieI].W() *edmFactor;
        }
        forAll(reactions_[reactionI].rhs() , rhsSpecieI)
        {
                label specieI = reactions_[reactionI].rhs()[rhsSpecieI].index;
                scalar stoichCoeff = reactions_[reactionI].rhs()[rhsSpecieI].stoichCoeff;
                RR_[specieI]+= stoichCoeff  * specieThermo_[specieI].W() *edmFactor;
        }
    }


    return 0;
}


// ************************************************************************* //
