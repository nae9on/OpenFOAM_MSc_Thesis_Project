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

#include "edmKinetic.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edmKinetic<CompType, ThermoType>::edmKinetic
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

    Info<< "edmKinetic: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edmKinetic<CompType, ThermoType>::~edmKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::scalarField Foam::edmKinetic<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalarField om((nSpecie_+2), 0.0);

    forAll(reactions_, i)
    {
        const Reaction<ThermoType>& R = reactions_[i];

        scalar omegai = omega
        (
            R, c, T, p, pf, cf, lRef, pr, cr, rRef
        );

        forAll(R.lhs(), s)
        {
            label si = R.lhs()[s].index;
            scalar sl = R.lhs()[s].stoichCoeff;
            om[si] -= sl*omegai;
        }

        forAll(R.rhs(), s)
        {
            label si = R.rhs()[s].index;
            scalar sr = R.rhs()[s].stoichCoeff;
            om[si] += sr*omegai;
        }
    }

    return om;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::edmKinetic<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c2(nSpecie_, 0.0);
    for (label i=0; i<nSpecie_; i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    scalar kf = R.kf(T, p, c2);
    scalar kr = R.kr(kf, T, p, c2);

    pf = 1.0;
    pr = 1.0;

    label Nl = R.lhs().size();
    label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        scalar exp = R.lhs()[slRef].exponent;
        if (exp<1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        scalar exp = R.rhs()[srRef].exponent;
        if (exp<1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edmKinetic<CompType, ThermoType>::tc() const
{
    FatalError << "edmKinetic::tc called, but not implemented!!"  << exit(FatalError);
    
    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edmKinetic<CompType, ThermoType>::dQ() const
{
    FatalError << "edmKinetic::dQ called, but not implemented!!" << exit(FatalError);

    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edmKinetic<CompType, ThermoType>::Sh() const
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
Foam::scalar Foam::edmKinetic<CompType, ThermoType>::solve
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

    // scalarField with laminar rates for one reaction
    PtrList<scalarField> RRLaminar(nSpecie_);
    forAll(Y_,specieI)
    {
        RRLaminar.set(specieI,new scalarField(rho.size(),0.));
    
    }    

    // scalarField with edm rates for one reaction  
    PtrList<scalarField> RREdm(nSpecie_);
    forAll(Y_,specieI)
    {
        RREdm.set(specieI,new scalarField(rho.size(),0.));
    }    

    // set all reaction rates to zero
    forAll(Y_,specieI)
    {
        RR_[specieI]=0.;
    }

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
                RREdm[specieI]-= stoichCoeff  * specieThermo_[specieI].W() *edmFactor;
        }
        forAll(reactions_[reactionI].rhs() , rhsSpecieI)
        {
                label specieI = reactions_[reactionI].rhs()[rhsSpecieI].index;
                scalar stoichCoeff = reactions_[reactionI].rhs()[rhsSpecieI].stoichCoeff;
                RREdm[specieI]+= stoichCoeff  * specieThermo_[specieI].W() *edmFactor;
        }
    

        //calculate laminar rates
    
        for(label celli=0; celli<rho.size(); celli++)
        {
            forAll(Y_,specieI)
            {
                RRLaminar[specieI][celli] = 0.0;
            }
            
            scalar rhoi = rho[celli];
            scalar Ti = this->thermo().T()[celli];
            scalar pi = this->thermo().p()[celli];
            scalarField c(Y_.size());
            scalarField dcdt(reactions_.size(), 0.0);
            forAll(Y_,specieI)
            {
                scalar Yi = Y_[specieI][celli];
                c[specieI] = rhoi*Yi/specieThermo_[specieI].W();
            }
            
            dcdt = omega(c, Ti, pi);
          
            forAll(Y_,specieI)
            {
                RRLaminar[specieI][celli] = dcdt[specieI]*specieThermo_[specieI].W();
            }
        }

        //combine laminar and edm rates, i.e. calculate harmonic mean
        forAll(Y_, specieI)
        {
            // prevent division by zero error if any reaction rate equals zero
            forAll(rho.internalField(), cellI)
            {
                if (RREdm[specieI][cellI]==0.)
                    RR_[specieI][cellI]=0.;
                else if (RRLaminar[specieI][cellI]==0.)
                    RR_[specieI][cellI]=0.;
                else 
                    // fill internal field with harmonic mean of chemical and EDM reaction rate         
                    RR_[specieI][cellI]+=
                               (RREdm[specieI][cellI] * RRLaminar[specieI][cellI])
                             / (RREdm[specieI][cellI] + RRLaminar[specieI][cellI] );

            }
           
        }        
    }
    
    return 0;
}


// ************************************************************************* //
