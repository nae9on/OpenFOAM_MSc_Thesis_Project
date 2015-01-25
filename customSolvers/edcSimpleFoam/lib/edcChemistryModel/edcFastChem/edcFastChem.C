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

#include "edcFastChem.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edcFastChem<CompType, ThermoType>::edcFastChem
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


    RR_(nSpecie_),
    YStar_(nSpecie_)
{
    // create the fields for the chemistry sources
    forAll(RR_, specieI)
    {
        RR_.set(specieI, new scalarField(mesh.nCells(), 0.));
        // allocate memory for fine structure data and initialize with cell average mass fractions
        YStar_.set(specieI, new scalarField( Y_[specieI].internalField() ));
        
    }

    // get access to chemistryPropertiesDict 
    const IOdictionary& chemistryProperties = this->db().objectRegistry::lookupObject<IOdictionary>("chemistryProperties");

    relaxFineStructures_=readScalar(chemistryProperties.lookup("underRelaxFineStructures"));

    // look for the names of the oxidiser specie in the subdict
    word oxidiserName=word(chemistryProperties.subDict("edcFastChemCoeffs").lookup("oxidiserName"));
    indexOxidiser_= dynamic_cast<const reactingMixture<ThermoType>&> (this->thermo()).species()[oxidiserName];

    // and for main fuel specie as well
    word mainFuelName=word(chemistryProperties.subDict("edcFastChemCoeffs").lookup("mainFuelName"));
    indexMainFuel_= dynamic_cast<const reactingMixture<ThermoType>&> (this->thermo()).species()[mainFuelName];

    // is additional output required?
    additionalOutput_=readBool(chemistryProperties.lookup("additionalOutput"));
    
    // find reaction which consumes main fuel
    forAll(reactions_, reactionI)
    {
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        { 
            label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;
             if (specieI==indexMainFuel_)
                 indexMainReaction_=reactionI;
        }
    }    

    Info << "edcFastChem: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edcFastChem<CompType, ThermoType>::~edcFastChem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcFastChem<CompType, ThermoType>::tc() const
{
    FatalError << "edcFastChem::tc called, but not implemented!!"  << exit(FatalError);
    
    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcFastChem<CompType, ThermoType>::dQ() const
{
    FatalError << "edcFastChem::dQ called, but not implemented!!" << exit(FatalError);

    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcFastChem<CompType, ThermoType>::Sh() const
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
void  Foam::edcFastChem<CompType, ThermoType>::writeScalarField(scalarField X, word fieldName)
{
    //output only if it is write time
    if ( this->mesh().time().outputTime() && additionalOutput_)
    {  
        // as paraView only undestands volScalarFields, make a volScalarField out of the scalarField and write it
        volScalarField volX
        (

                IOobject
                (
                    fieldName,
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
        );

        volX.internalField() = X;
        volX.correctBoundaryConditions();
        volX.write();    
    }
}


template<class CompType, class ThermoType>
Foam::PtrList<Foam::scalarField >  Foam::edcFastChem<CompType, ThermoType>::calculateYStar()
{
    // normalized mass fractions for resulting reactor state
    PtrList<scalarField> YStar(Y_.size());
    forAll(YStar, specieI)
    {
        // construct YStar fields as copies from Y_
        YStar.set(specieI, new scalarField( Y_[specieI].internalField()  ) );
    }

    // normalization factor = stoichCoeff * molecularWeight
    List<scalar> normFactor(Y_.size());

    forAll(reactions_ , reactionI)
    {
    //Unsolved: What happenes if a species is on both sides of the reaction ???
    //But this should not be the case in a realistic global chemical reaction.

        normFactor = 0.;
        // calculate norm factor for reaction with index: reactionI
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        {
                label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;
                scalar stoichCoeff = reactions_[reactionI].lhs()[lhsSpecieI].stoichCoeff;
                normFactor[specieI]-=stoichCoeff  * specieThermo_[specieI].W();            
        }
        forAll(reactions_[reactionI].rhs() , rhsSpecieI)
        {
                label specieI = reactions_[reactionI].rhs()[rhsSpecieI].index;
                scalar stoichCoeff = reactions_[reactionI].rhs()[rhsSpecieI].stoichCoeff;
                normFactor[specieI]+=stoichCoeff  * specieThermo_[specieI].W();            
        }

        // normalize all mass fractions according to abs(normFactor)
        forAll(Y_, specieI)
        {

            normFactor[specieI]=fabs(normFactor[specieI]);
            // if stoichCoeff ist zero => normFactor is zero: in this case prevent division by zero
            // do not normalize if normFactor is zero (in this case YStar will be untouched until the end the loop iteration)
            if (normFactor[specieI] > SMALL)
                YStar[specieI]/=normFactor[specieI];
            Info << "normFactor(" << specieI<<") = "<<normFactor[specieI];
        }
        

        // add all fuel mass fractions (i.e. all on LHS exept for oxidiser)

        // normalized mass fraction of all the fuels combined
        scalarField YStarFuel(YStar[0].size(),0.);
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        {            
            label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;    
            scalar stoichCoeff = reactions_[reactionI].lhs()[lhsSpecieI].stoichCoeff;
            // add to fuel only if the stoichCoeff is not zero (and of course it is not oxidiser)        
            if ((specieI!= indexOxidiser_) && ( stoichCoeff > SMALL))
            {
                YStarFuel+=YStar[specieI];
            }
        }

        // normalized mass fraction with oxidiser
        scalarField YStarOxidiser=YStar[indexOxidiser_];

        // calculate reactor state for fast chemistry
        //how much can react at all?
        scalarField reactedMass=min(YStarFuel,YStarOxidiser);

        // comsume reactants
        // as all fuel has been sumarized to YStarFuel, consume is weighted by mass fraction
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        {
                label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;
                scalar stoichCoeff = reactions_[reactionI].lhs()[lhsSpecieI].stoichCoeff;
                if ((specieI!= indexOxidiser_) && ( stoichCoeff != 0))  // this is a fuel specie
                {
                    YStar[specieI]-=(reactedMass * YStar[specieI]/max(YStarFuel,SMALL));
                }
                else if (specieI==indexOxidiser_)// this is the oxidiser
                {
                    YStar[indexOxidiser_]-=reactedMass;
                }
        }
        //build up products
        forAll(reactions_[reactionI].rhs() , rhsSpecieI)
        {
                label specieI = reactions_[reactionI].rhs()[rhsSpecieI].index;
                YStar[specieI]+=reactedMass;
        }         
        
        // de-normalize mass fractions
        forAll(YStar, specieI)
        {
            if (normFactor[specieI] > SMALL)
                YStar[specieI]*=normFactor[specieI];
        }        
        
        // and do this again for the next reaction (which has different normalization factors)        
    }
    return YStar;
}

template<class CompType, class ThermoType>
Foam::scalar Foam::edcFastChem<CompType, ThermoType>::solve
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
    Info << "Selected chemistry model is EDC/FastChemistry" << endl;

    // set all reaction rates to zero
    forAll(Y_,specieI)
    {
        RR_[specieI]=0.;
    }    
   
    // get flow parameters: rho, k, epsilon, nu    
    scalarField& rho = this->thermo().rho().internalField();
    const scalarField& k = this->db().objectRegistry::lookupObject<volScalarField>("k").internalField();         
    scalarField epsilon(rho.size(),0.);
    if (this->db().objectRegistry::foundObject<volScalarField>("epsilon"))
    {
	const scalarField& epsilonConst = this->db().objectRegistry::lookupObject<volScalarField>("epsilon").internalField();
	epsilon = epsilonConst;
    }
    // or maybe we can get omega 
    else if (this->db().objectRegistry::foundObject<volScalarField>("omega"))
    {
	const scalarField& omega = this->db().objectRegistry::lookupObject<volScalarField>("omega").internalField();
        //calculate epsilon value
	const scalar Cmu=0.09;
	epsilon = Cmu*omega*k;

    }    
    // or we are lost!!
    else
    {
	FatalError << "Could not access epsilon or omega data from turbulence model!!" << exit(FatalError);
    }                     
    scalarField nu=this->thermo().mu().internalField()/rho;

    
    // calculate EDC model variables: gammaStar and mDotStar
    scalarField gammaL= 2.1377 * pow(  nu*epsilon / pow(k,2.)  , (1./4.) );
    gammaL= min (0.7, gammaL);
    scalarField mDotStar= 2.45*pow( epsilon/nu  , (1./2.) );

    // calculate fine structure composition according to fast chemistry approach
    PtrList<scalarField> YStarNew=calculateYStar();
    // do explicit underrelaxation of fine structures
    forAll(YStar_,specieI)
    {
        YStar_[specieI]=YStar_[specieI]*(1.-relaxFineStructures_) + YStarNew[specieI]*relaxFineStructures_;
    }
    

    //PtrList to normalized mass fractions
    PtrList<scalarField> YNorm(Y_.size());
    
    // normalization factor = stoichCoeff * molecularWeight
    List<scalar> normFactor(Y_.size());

    normFactor = 0.;
    // calculate norm factor for main fuel reaction
    forAll(reactions_[indexMainReaction_].lhs() , lhsSpecieI)
    {
            label specieI = reactions_[indexMainReaction_].lhs()[lhsSpecieI].index;
            scalar stoichCoeff = reactions_[indexMainReaction_].lhs()[lhsSpecieI].stoichCoeff;
            normFactor[specieI]-=stoichCoeff  * specieThermo_[specieI].W();            
    }
    forAll(reactions_[indexMainReaction_].rhs() , rhsSpecieI)
    {
            label specieI = reactions_[indexMainReaction_].rhs()[rhsSpecieI].index;
            scalar stoichCoeff = reactions_[indexMainReaction_].rhs()[rhsSpecieI].stoichCoeff;
            normFactor[specieI]+=stoichCoeff  * specieThermo_[specieI].W();            
    }
    
    // normalize mass fractions according to normFactor
    forAll(YStar_, specieI)
    {
        normFactor[specieI]= fabs(normFactor[specieI]);
        // construct YStar_ fields as copies from Y_
        YNorm.set(specieI, new scalarField( Y_[specieI].internalField() ));
        // and normalize them
        if (normFactor[specieI] > SMALL)
            YNorm[specieI] /= normFactor[specieI];
    }

    // create "fuel" mass fraction from all species on LHS,
    // exept for oxidiser and species with zero stoichiometric coefficient
    scalarField YNormFuel(rho.size(), 0.);
    forAll(reactions_[indexMainReaction_].lhs() , lhsSpecieI)
    {            
        label specieI = reactions_[indexMainReaction_].lhs()[lhsSpecieI].index;    
        scalar stoichCoeff = reactions_[indexMainReaction_].lhs()[lhsSpecieI].stoichCoeff;
        // add to fuel only if the stoichCoeff is not zero (and of course it is not oxidiser)        
        if ((specieI!= indexOxidiser_) && ( stoichCoeff > SMALL))
        {
            YNormFuel+=YNorm[specieI];
        }
    }        

    //  create "products" mass fraction from reaction equation RHS
    scalarField YProducts(rho.size(), 0.);
    forAll(reactions_[indexMainReaction_].rhs() , rhsSpecieI)
    {
        label specieI = reactions_[indexMainReaction_].rhs()[rhsSpecieI].index;
        scalar stoichCoeff = reactions_[indexMainReaction_].rhs()[rhsSpecieI].stoichCoeff;
        // if stoichCoeff is zero, don't consider the species' mass fraction
        if (stoichCoeff > SMALL)
            YProducts+=YNorm[specieI];
    }

    //  look for minimum value in reaction equations LHS
    scalarField YMin=min(YNorm[indexOxidiser_], YNormFuel);

    scalarField X1= pow(YMin+YProducts, 2.) / max(YNormFuel+YProducts,SMALL) / max(YNorm[indexOxidiser_]+YProducts,SMALL);

    scalarField X2=  1/gammaL * YProducts/ max(YProducts+YMin,SMALL);
    X2= min(1., X2);

    scalarField X3=  gammaL * (YProducts+YMin)/ max(YMin,SMALL) ;
    X3= min(1., X3);
            
    scalarField X= X1*X2*X3; 

    
    // calculate reaction rate coefficient                                          
    scalarField edcFactor = rho * pow(gammaL,3.) * X * mDotStar / (1-pow(gammaL,3.)*X);

    //output some fields (if requested)
    writeScalarField(gammaL,"gammaL");
    writeScalarField(X,"X");
    writeScalarField(mDotStar,"mDotStar");




                    
                    
    //  transform reaction rate to un-normalized specie space
    forAll(Y_,specieI)
    {
        RR_[specieI]= (-1) * edcFactor*(Y_[specieI]-YStar_[specieI]);  
    }
    
    
    //output fine structure data 
    forAll(Y_,specieI)
    {
        writeScalarField(YStar_[specieI],"YStar(" + this->Y_[specieI].name() + ")");
    }            
    
    return 0;
}


// ************************************************************************* //
