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

#include "edcLE.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edcLE<CompType, ThermoType>::edcLE
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
    YStar_(nSpecie_),
    LETemp_(0),
    tauChMin_(0),
    tauCh_(mesh.nCells(),0.)
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

    // read Local Extinction polynomial coefficients from chemistryProperties
    relaxFineStructures_=readScalar(chemistryProperties.lookup("underRelaxFineStructures"));

    // is additional output required?
    additionalOutput_=readBool(chemistryProperties.lookup("additionalOutput"));


    // read Local Extinction polynomial coefficients from chemistryProperties
    label curveI = 1;
    word curveLabel= "curve" + itoa(curveI) ;

    while ( chemistryProperties.subDict("edcLECoeffs").found(curveLabel) )
    {
        // resize lists
        LETemp_.setSize(curveI);
        tauChMin_.setSize(curveI);

        
        // read curve values from chemistryProperties and put them in the corresponding list
        LETemp_[curveI-1]=readScalar(chemistryProperties.subDict("edcLECoeffs").subDict(curveLabel).lookup("temperature"));
        tauChMin_[curveI-1]=readScalar(chemistryProperties.subDict("edcLECoeffs").subDict(curveLabel).lookup("tauChMin"));        
        
        curveI++;
        curveLabel= "curve" + itoa(curveI) ;
    }
  
    // and a litle error checking for empty lists
    if (curveI == 1)
    {
        FatalError << "No local extinction curves defined in edcLECoeffs!!"  << exit(FatalError);
    }
    
    // look for the names of the oxidiser specie in the subdict
    word oxidiserName=word(chemistryProperties.subDict("edcLECoeffs").lookup("oxidiserName"));
    indexOxidiser_= dynamic_cast<const reactingMixture<ThermoType>&> (this->thermo()).species()[oxidiserName];
    
    // and for main fuel specie as well
    word mainFuelName=word(chemistryProperties.subDict("edcLECoeffs").lookup("mainFuelName"));
    indexMainFuel_= dynamic_cast<const reactingMixture<ThermoType>&> (this->thermo()).species()[mainFuelName];

    //read auto-ignition temperature from chemistryProperties
    autoIgnitionTemperature_=readScalar(chemistryProperties.subDict("edcLECoeffs").lookup("autoIgnitionTemperature"));

    Info << "Species defined as oxidiser in control/chemistryPropeties: " << oxidiserName << endl << endl;
    Info << "Fuel consumption priorities: " << endl;     
    

    
    fStoich_=0.;
    forAll(reactions_, reactionI)
    {
        forAll(reactions_[reactionI].lhs() , lhsSpecieI)
        {            
            label specieI = reactions_[reactionI].lhs()[lhsSpecieI].index;
            scalar stoichCoeff = reactions_[reactionI].lhs()[lhsSpecieI].stoichCoeff;
            
            if ((specieI != indexOxidiser_) && (stoichCoeff > SMALL))
                Info << dynamic_cast<const reactingMixture<ThermoType>&> (this->thermo()).species()[specieI];
                
            // does this reaction contain the main fuel?            
            if ((specieI==indexMainFuel_) && (fStoich_==0.))
            {
                scalar stoichMassFuel=0.;
                scalar stoichMassOxidiser=0.;
                
                forAll(reactions_[reactionI].lhs() , thislhsSpecieI)
                {            
                    label thisspecieI = reactions_[reactionI].lhs()[thislhsSpecieI].index;
                    scalar thisstoichCoeff = reactions_[reactionI].lhs()[thislhsSpecieI].stoichCoeff;
                    if ((thisspecieI != indexOxidiser_) && (thisstoichCoeff > SMALL)) // is this a fuel component?
                    {
                        //yes, so let's add it up ( += stoichCoeff * molarWeight)
                        stoichMassFuel += thisstoichCoeff * specieThermo_[thislhsSpecieI].W();
                    }
                    if ((thisspecieI == indexOxidiser_) && (thisstoichCoeff > SMALL)) // is this oxidiser?
                    {
                        //yes, so let's add it up ( += stoichCoeff * molarWeight)
                        stoichMassOxidiser += thisstoichCoeff * specieThermo_[thislhsSpecieI].W();
                    }
                    
                }
  
                // stochiometric mass fuel/ mass air ratio 
                fStoich_ = stoichMassFuel / stoichMassOxidiser;                
                
                indexMainReaction_=reactionI;
            }
        }
    

        Info << endl;
    }
    Info << "Stoiciometric fuel ratio calculated from reaction #" << indexMainReaction_ << ": " << fStoich_ << endl;
    Info << endl;
    Info << "edcLE: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edcLE<CompType, ThermoType>::~edcLE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcLE<CompType, ThermoType>::tc() const
{
    FatalError << "edcLE::tc called, but not implemented!!"  << exit(FatalError);
    
    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcLE<CompType, ThermoType>::dQ() const
{
    FatalError << "edcLE::dQ called, but not implemented!!" << exit(FatalError);

    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcLE<CompType, ThermoType>::Sh() const
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
bool  Foam::edcLE<CompType, ThermoType>::isLocalExtinction(Foam::scalar temp, Foam::scalar equivalenceRatio, Foam::scalar tauEdc)
{
    reason_=0.;
    if( LETemp_.size() == 1)
    {
        // check tauCh > tauEdc (=> local extinciton occurs)
        if ( tauChMin_[0] > tauEdc)
            return true;
    }
    else if ( LETemp_.size() > 1)
    {
        // if we are above autoignition temperature, there is no local extinction
        if (temp > autoIgnitionTemperature_)
        {
            return false;
        }
        else
        {
            //if not, check for tauCh 
    
            scalar tauChMin=interpolateXY(temp,LETemp_,tauChMin_);
            // check tauCh > tauEdc (=> local extinciton occurs)
            if ( tauChMin > tauEdc)
            {
                reason_=-1.;
                return true;
            }        
        
        }
                
            
    }
    else
        FatalError << "Error in edcLE::isLocalExtinction(): LE data lists are emtpy!"  << exit(FatalError);

    // otherwise no extinction
    return false;
        
    
}

template<class CompType, class ThermoType>
void  Foam::edcLE<CompType, ThermoType>::writeScalarField(scalarField X, word fieldName)
{
    //output only if it is write time and "additionalOutput" is requested
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
Foam::PtrList<Foam::scalarField >  Foam::edcLE<CompType, ThermoType>::calculateYStar()
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
    //Unresolved: What happenes if a species is on both sides of the reaction ???
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
Foam::scalar Foam::edcLE<CompType, ThermoType>::solve
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
    Info << "Selected chemistry model is EDC/LocalExtinction" << endl;

    // set all reaction rates to zero
    forAll(Y_,specieI)
    {
        RR_[specieI]=0.;
    }    
   
    // get flow parameters: rho, k, epsilon, nu    
    scalarField& rho = this->thermo().rho().internalField();
    const scalarField& temperature = this->thermo().T().internalField();
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
    // max(k,0.) to prevent FPE is k < 0 (which may happen if k-epsilon has problems)
    scalarField gammaL= 2.1377 * pow(  nu*epsilon / pow(max(k,0.),2.)  , (1./4.) );
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

    writeScalarField(gammaL,"gammaL");
    writeScalarField(X,"X");
    writeScalarField(mDotStar,"mDotStar");



    //Calculate LOCAL EXTINCTION (for main fuel reaction only)

    // calculate equivalence ratio
    scalarField equivalenceRatio;
    equivalenceRatio = Y_[indexMainFuel_] / max(Y_[indexOxidiser_],SMALL) / fStoich_;
    writeScalarField(equivalenceRatio,"eR");


    scalarField localExtinction(rho.size(),0.);
    
    forAll(rho,cellI)
    {
         scalar T=temperature[cellI];
         scalar eR=equivalenceRatio[cellI];

        if ( isLocalExtinction(T,eR,1/mDotStar[cellI]))
        {
            edcFactor[cellI]=0.;
            localExtinction[cellI]=reason_;
        }
        else
             localExtinction[cellI]=1;           
             
   }

    writeScalarField(localExtinction,"localExt");
                    
                    
    //  transform reaction rate to un-normalized specie space
    forAll(Y_,specieI)
    {
        RR_[specieI]= (-1) * edcFactor*(Y_[specieI]-YStar_[specieI]);  
    }
    
    //output fine structure data 
    if ( this->mesh().time().outputTime())
    {         
        forAll(Y_,specieI)
        {
            writeScalarField(YStar_[specieI],"YStar(" + this->Y_[specieI].name() + ")");
        }            
    } 

    return 0;
}


// ************************************************************************* //
