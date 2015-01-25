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

#include "edcPSR.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edcPSR<CompType, ThermoType>::edcPSR
(
    const fvMesh& mesh,
    const word& compTypeName,
    const word& thermoTypeName
)
:
    CompType(mesh, thermoTypeName),
    // initialize RADAU5 parameters:
	StiffIntegratorT(
	    this->thermo().composition().Y().size()+1,      //number of equations to solve (=number of species + 1)

	    solverX,            //solution vector (double array of dimension nSpecie_+1)
	    0.,                 //startTime (will be reset later)
	    10.,                //endTime (will be reset later)
		1.,                 //time step for intermediate output
		0,                  //itoler = 0:	Both rtoler and atoler are scalars; itoler = 1:	Both rtoler and atoler are array
	    new scalar(readScalar((this->db().objectRegistry::lookupObject<IOdictionary>("chemistryProperties")).subDict("edcPSRCoeffs").lookup("relativeTolerance"))),   
		                    //relative tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
		new scalar(readScalar((this->db().objectRegistry::lookupObject<IOdictionary>("chemistryProperties")).subDict("edcPSRCoeffs").lookup("absoluteTolerance"))),
		                    //absolute tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
		false,               //provide data for output functions
		1.e-3,                 //initial step size guess (usually 1.0e-3 or 1.0e-5; not very important!; if 0. code sets h = 1.0e-6)
		0.,                 //maximal step size ( if 0. code sets (xend - x) )
		readScalar((this->db().objectRegistry::lookupObject<IOdictionary>("chemistryProperties")).subDict("edcPSRCoeffs").lookup("maxIterations")) ,         
		                    //max. number of steps ( if 0. then 1e5 is assumed)
		0.,                 //rounding unit (if 0. then 1.0e-16 is assumed)
		0.,                 //The safety factor in step size prediction ( if 0. then 0.9 is assumed)
		0.,                 //Parameter for step size selection (if 0. then fac1 = 5.0 is assumed), 
		0.,                 //Parameter for step size selection (if 0. then facr = 1.0/8.0 is assumed)
		false,              //Switch for the computation of the Jacobian (0= Jacobian calculated internally)
		this->thermo().composition().Y().size()+1,
		                    //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
		this->thermo().composition().Y().size()+1,        
		                    //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
		0,                  //Gives information on the mass-matrix ( 0 = mass-matrix is unity)
		0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
		0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
		0,                  //maximal number of Newton iterations ( if 0 then 7 is assumed)
		true,              //If startn != 0 zero starting values are used (beneficial if convergence problems occur)
		0,                  //Dimension of the index 1 variables (must be > 0; if 0 then nSpecie_+1 is assumed)
		0,                  //Dimension of the index 2 variables. Default nind2 = 0.
		0,                  //Dimension of the index 2 variables. Default nind2 = 0.
		1,                  //Switch for step size strategy: npred = 1 or 0  mod. predictive controller (Gustafsson) [default]; 
		                    //                            If npred = 2  classical step size control (less safe; slightly faster)
		0,                  //Default m1 = 0
		0,                  //Default m2 = m1
		false,              //If != 0 Jacobian matrix transformed to Hessenberg form 
		                    //(advantageous for large systems with full Jacobian)
		0,                  //Stopping criterion for Newton's method, usually chosen < 1. 
		                    //Smaller values of fnewt make the code slower, but safer; Default min(0.03, sqrt(rtoler))
		0.,                 //If quot1 < hnew/hold < quot2, then the step size is not changed.
		0.,                 //This saves, together with a large thet, lu-decompositions and
            				//computing time for large systems. for small systems one may have
            				//quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
            				//quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

		0.),                //Criterion for recomputing Jacobian; default 0.001; Increase thet to 0.1 if Jacobian costly


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
    mDotStar_(mesh.nCells(),0.),
    gammaStar_(mesh.nCells(),0.),
    solverX(new double[nSpecie_+1])

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

    // ODE integration absolute tolerance
    absoluteTolerance_=readScalar(chemistryProperties.subDict("edcPSRCoeffs").lookup("absoluteTolerance"));    

    // ODE integration relative tolerance
    relativeTolerance_=readScalar(chemistryProperties.subDict("edcPSRCoeffs").lookup("relativeTolerance"));    

    // ODE integration max. allowed number of iterations
    maxIterations_=readScalar(chemistryProperties.subDict("edcPSRCoeffs").lookup("maxIterations"));

    // use binaryTable to store and retrieve ODE integration results?
    useBinaryTree_=readBool(chemistryProperties.subDict("edcPSRCoeffs").lookup("useBinaryTree"));    

    // tolerance for retrieving values and growing points
    tableErr_=readScalar(chemistryProperties.subDict("edcPSRCoeffs").lookup("binaryTreeTolerance"));    

    // tolerance for retrieving values and growing points (read as scalar to allew e.g. "1.e5" in chemistryDict)
    tableSize_=readScalar(chemistryProperties.subDict("edcPSRCoeffs").lookup("binaryTreeSize"));    
    
    
    TMin_=dimensionedScalar(this->thermo().lookup("TMin")).value();
    TMax_=dimensionedScalar(this->thermo().lookup("TMax")).value();    

    
    // initialize binaryTree (if wanted)
    if(useBinaryTree_)
    {
        resultTable_=new binaryTree(label(tableSize_));

        
        Info << "maxTableSize: " << label(tableSize_) << endl;

        //set these values to <0 to indicate that they have not yet been initialized
        hsMax_ = -1;
        mDotStarMax_ = -1;
        gammaStarMax_ = -1;        
    }
    

    Info << "edcPSR: Number of species = " << nSpecie_
        << " and reactions = " << nReaction_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::edcPSR<CompType, ThermoType>::~edcPSR()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::scalarField Foam::edcPSR<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    scalarField om((nSpecie_), 0.0);

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
Foam::scalar Foam::edcPSR<CompType, ThermoType>::omega
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
Foam::edcPSR<CompType, ThermoType>::tc() const
{
    FatalError << "edcPSR::tc called, but not implemented!!"  << exit(FatalError);
    
    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}


template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcPSR<CompType, ThermoType>::dQ() const
{
    FatalError << "edcPSR::dQ called, but not implemented!!" << exit(FatalError);

    //the return statement should never get called
    volScalarField dummyField(Y_[0]);
    return dummyField;
}

template<class CompType, class ThermoType>
Foam::tmp<Foam::volScalarField>
Foam::edcPSR<CompType, ThermoType>::Sh() const
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
void  Foam::edcPSR<CompType, ThermoType>::writeScalarField(scalarField X, word fieldName)
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
void Foam::edcPSR<CompType, ThermoType>::iterateTAndRho(Foam::scalar& TResult, Foam::scalar& rho, Foam::scalar h, Foam::scalar p, Foam::scalarField& Y, Foam::scalar TInit) const
{
    ThermoType mixture(0.0*specieThermo_[0]);
    // use size from global variable, so that YStarCell (which contains the cell index as Y_.size+1) 
    forAll(Y_,specieI)
    {
	    mixture += Y[specieI]/specieThermo_[specieI].W()*specieThermo_[specieI];
    }     
    
     TResult = TInit;
     scalar Ttol = TInit*1e-4;
     label iter = 0;

     //Newton iteration
     do
     {
        
         TInit = TResult;
         TResult = TInit - (mixture.Hs(TInit) - h)/mixture.Cp(TInit);
         //clip value since calling H(T) and Cp(T) out of range is a critical error
         TResult=max(TResult,TMin_);
	     TResult=min(TResult,TMax_);

    	 iter++;
 
     } while ((mag(TResult - TInit) > Ttol) && (iter < 100));


     if (iter >= 100)
     {
	    TResult = TInit;
	    Info << "Max number of iterations reached; setting TResult to " << TInit << endl;
     }


    // calculate mean mole fraction in fine structures
    scalar MMean=0.;
    forAll(Y_,specieI)
    {
	    MMean += Y[specieI]/specieThermo_[specieI].W();
    }
    MMean = 1./MMean;

    rho= p*MMean/ (TResult*mixture.RR);
    
}

template<class CompType, class ThermoType>
void Foam::edcPSR<CompType, ThermoType>::derivative
(Foam::scalarField&  YStarCell, Foam::scalarField& dYdtCell) const
{    
    //which cell are we solving for?
    scalar cellI = YStarCell[nSpecie_];

    // surrounding fluid mass fractions
    scalarField YSurrCell(nSpecie_);


    forAll(YStar_,specieI)
    {
        // clip YStarCell to positive values
        YStarCell[specieI]=max(YStarCell[specieI],0.0);
        // calculate surrounding's mass fractions
        YSurrCell[specieI]=(Y_[specieI].internalField()[cellI] - YStarCell[specieI] * gammaStar_[cellI]) 
                                        / ( 1.- gammaStar_[cellI]);
    }

    // assume hStar = hBar !!
    scalar hStarCell=this->thermo().hs().internalField()[cellI];

    // pStar = pBar !!
    scalar pCell=this->thermo().p().internalField()[cellI];

    //update temperature    
    scalar TStarCell = 0.;
    scalar rhoStarCell = 0.;
    iterateTAndRho(TStarCell, rhoStarCell, hStarCell, pCell , YStarCell, this->thermo().T().internalField()[cellI]);
    
    // now calculate fine structure concentrations with rhoStar
    scalarField CStarCell(nSpecie_);
    //initialize YStarCell and YSurrCell from y
    forAll(YStar_,specieI)
    {
        CStarCell[specieI]=max(YStarCell[specieI]*rhoStarCell/specieThermo_[specieI].W(),1e-20);
    }   
    //calculate chemistry kinetics    

    scalarField dCdtCell(nSpecie_);
    
    dCdtCell = omega(CStarCell,TStarCell,pCell);


    //convert concentrations to mass fractions
    // dYdt is chemical kinetics plus m* (Y°-Y*)
    forAll(YStar_,specieI)
    {
        dYdtCell[specieI]=dCdtCell[specieI]/rhoStarCell*specieThermo_[specieI].W() 
                                    +  mDotStar_[cellI] * (YSurrCell[specieI]-YStarCell[specieI]) ;
    }  
    
    
    dYdtCell[nSpecie_]= 0.; // we are not solving for the cell index



}


template<class CompType, class ThermoType>
void  Foam::edcPSR<CompType, ThermoType>::updateYStar(scalar startTime, scalar stopTime)
{
    if(useBinaryTree_)
    {
        //USE BINARY TABULATION
        
        //initialize some statistical data
        label added=0;
        label grown=0;
        label retrieved=0;
       
        // if not yet happened, initialize normalization values for hs, mDotStar and gammaStar
        // (NB this cannot be done in constructor, since there is no turbulence yet defined [requires thermo first])

        if(mDotStarMax_ < 0)
        {
            //get max order of magnitude of some fields
            hsMax_ = Foam::max(mag(this->thermo().hs())).value();
            mDotStarMax_ = Foam::max(mDotStar_);
            gammaStarMax_ = Foam::max(gammaStar_);  
        }   
        
        forAll(Y_[0].internalField(), cellI)
        {       
            // initialize query vector (Yi , hStar/hBar , mDotStar, gammStar)
            scalarField queryVector(Y_.size()+3);
            //load speciesData
            forAll(Y_,specieI)
            {
                queryVector[specieI]=YStar_[specieI][cellI];        
            }
            
            // Info << " max Values" << hsMax_ << "   " << mDotStarMax_ << "   " << gammaStarMax_ << endl;
            
            // add further values and normalize them to more or less 0..1 
            queryVector[Y_.size()] = this->thermo().hs().internalField()[cellI]/hsMax_;            
            queryVector[Y_.size()+1] = mDotStar_[cellI]/mDotStarMax_;
            queryVector[Y_.size()+2] = gammaStar_[cellI]/gammaStarMax_;                           
            
            //query reactionTable
            chemPoint* tableEntry = resultTable_->findClosest(queryVector);

            //found something??
            if (tableEntry != NULL && ( tableEntry->inEOA(queryVector) ))
            {

                scalarField resultVector = tableEntry->r();
                //unpack result vector
                //unpack species (and underrelax)
                forAll(YStar_,specieI)
                {
                    YStar_[specieI][cellI]=YStar_[specieI][cellI]*(1.-relaxFineStructures_) + resultVector[specieI]*relaxFineStructures_;

                }
                
                // some statistical info
                retrieved++;
                //and continue with next cell
                continue;
            }
            else
            {
                // retrieve unsuccessfull, do direct integration with and add result to table (or grow it)
                //initialize solution vector
                forAll(Y_,specieI)
                {
                    solverX[specieI]=YStar_[specieI][cellI];        
                }
                solverX[nSpecie_]=cellI;

                // set stopTime for this cell
                stopTime=min((stopTime-startTime),100*1/(max(mDotStar_[cellI],SMALL)));

                // for this, reset the RADAU5 parameters
                double atol = absoluteTolerance_;
                double rtol = relativeTolerance_;

                // tell RADAU5 the start and endTimes
            	ResetT(
            	    solverX,            //solution vector (double array of dimension nSpecie_+1)
            	    startTime,                 //startTime (will be reset later)
            	    stopTime,                //endTime (will be reset later)
            		1.,                 //time step for intermediate output
            		0,                  //itoler = 0:	Both rtoler and atoler are scalars; itoler = 1:	Both rtoler and atoler are array
            	    &atol,   
            		                    //relative tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
            		&rtol,
            		                    //absolute tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
            		0.,                 //initial step size guess (usually 1.0e-3 or 1.0e-5; not very important!; if 0. code sets h = 1.0e-6)
            		0.,                 //maximal step size ( if 0. code sets (xend - x) )
            		maxIterations_ ,         
            		                    //max. number of steps ( if 0. then 1e5 is assumed)
            		0.,                 //rounding unit (if 0. then 1.0e-16 is assumed)
            		0.,                 //The safety factor in step size prediction ( if 0. then 0.9 is assumed)
            		0.,                 //Parameter for step size selection (if 0. then fac1 = 5.0 is assumed), 
            		0.,                 //Parameter for step size selection (if 0. then facr = 1.0/8.0 is assumed)
            		nSpecie_+1,
            		                    //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
            		nSpecie_+1,        
            		                    //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)

            		0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
            		0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
            		0,                  //maximal number of Newton iterations ( if 0 then 7 is assumed)
            		true,              //If startn != 0 zero starting values are used (beneficial if convergence problems occur)
            		0,                  //Dimension of the index 1 variables (must be > 0; if 0 then nSpecie_+1 is assumed)
            		0,                  //Dimension of the index 2 variables. Default nind2 = 0.
            		0,                  //Dimension of the index 2 variables. Default nind2 = 0.
            		1,                  //Switch for step size strategy: npred = 1 or 0  mod. predictive controller (Gustafsson) [default]; 
            		                    //                            If npred = 2  classical step size control (less safe; slightly faster)
            		0,                  //Default m1 = 0
            		0,                  //Default m2 = m1
            		false,              //If != 0 Jacobian matrix transformed to Hessenberg form 
            		                    //(advantageous for large systems with full Jacobian)
            		0,                  //Stopping criterion for Newton's method, usually chosen < 1. 
            		                    //Smaller values of fnewt make the code slower, but safer; Default min(0.03, sqrt(rtoler))
            		0.,                 //If quot1 < hnew/hold < quot2, then the step size is not changed.
            		0.,                 //This saves, together with a large thet, lu-decompositions and
                        				//computing time for large systems. for small systems one may have
                        				//quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
                        				//quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

            		0.),                //Criterion for recomputing Jacobian; default 0.001; Increase thet to 0.1 if Jacobian costly
                // and solve ODEs for this cell

                Integrate();        
                // => solverX now contains result of ODE integration

                // do explicit underrelaxation of fine structures
                forAll(YStar_,specieI)
                {
                    YStar_[specieI][cellI]=YStar_[specieI][cellI]*(1.-relaxFineStructures_) + solverX[specieI]*relaxFineStructures_;
                }
                
                //store result in table
                scalarField resultVector(Y_.size());
                forAll(Y_,specieI)
                {
                    resultVector[specieI]=YStar_[specieI][cellI];        
                }

                scalarField  queryTolerances(queryVector.size(),1.);
                scalarField  resultTolerances(resultVector.size(),1.);            

                // check if grow is possible
                if (tableEntry!=NULL)
                {                
                    if (!tableEntry->checkSolution(queryVector,resultVector))
                    {
                        //not possible, add new leaf
                        chemPoint newEntry(queryVector,resultVector,queryTolerances,resultTolerances,tableErr_);
                        resultTable_->insert(newEntry);
                        added ++;
                    } 
                    else
                    {
                        grown++;
                    }
                }
                else
                {
                    chemPoint newEntry(queryVector,resultVector,queryTolerances,resultTolerances,tableErr_);
                    resultTable_->insert(newEntry);
                    added ++;
                }
            }   // end of integration
                        
        }   // end of cell loop

        //and output some statistics
        Info << "table size: " << resultTable_->size() << endl;
        Info << "added points: " << added << endl;
        Info << "retrieved points: " << retrieved << endl;    
        Info << "grown points: " << grown << endl;    

       
    }
    else
    {
        //DO NOT USE TABULATION, INTEGRATE FOR EVERY CELL
        forAll(Y_[0].internalField(), cellI)
        {
            //initialize solution vector
            forAll(Y_,specieI)
            {
                solverX[specieI]=YStar_[specieI][cellI];        
            }
            solverX[nSpecie_]=cellI;
            

            // set stopTime for this cell
            stopTime=min((stopTime-startTime),100*1/(max(mDotStar_[cellI],SMALL)));

            // for this, reset the RADAU5 parameters
            double atol = absoluteTolerance_;
            double rtol = relativeTolerance_;

            // tell RADAU5 the start and endTimes
        	ResetT(
        	    solverX,            //solution vector (double array of dimension nSpecie_+1)
        	    startTime,                 //startTime (will be reset later)
        	    stopTime,                //endTime (will be reset later)
        		1.,                 //time step for intermediate output
        		0,                  //itoler = 0:	Both rtoler and atoler are scalars; itoler = 1:	Both rtoler and atoler are array
        	    &atol,   
        		                    //relative tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
        		&rtol,
        		                    //absolute tolerance (either one value or array of size nSpecie_+1; if 0. then 1e-7 is assumed)
        		0.,                 //initial step size guess (usually 1.0e-3 or 1.0e-5; not very important!; if 0. code sets h = 1.0e-6)
        		0.,                 //maximal step size ( if 0. code sets (xend - x) )
        		maxIterations_ ,         
        		                    //max. number of steps ( if 0. then 1e5 is assumed)
        		0.,                 //rounding unit (if 0. then 1.0e-16 is assumed)
        		0.,                 //The safety factor in step size prediction ( if 0. then 0.9 is assumed)
        		0.,                 //Parameter for step size selection (if 0. then fac1 = 5.0 is assumed), 
        		0.,                 //Parameter for step size selection (if 0. then facr = 1.0/8.0 is assumed)
        		nSpecie_+1,
        		                    //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)
        		nSpecie_+1,        
        		                    //Switch for the banded structure of the Jacobian ( for details see RADAU5/StiffIntegratorT.cpp)

        		0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
        		0.,                 //Switch for the banded structure of the mass-matrix ( for details see RADAU5/StiffIntegratorT.cpp)
        		0,                  //maximal number of Newton iterations ( if 0 then 7 is assumed)
        		true,              //If startn != 0 zero starting values are used (beneficial if convergence problems occur)
        		0,                  //Dimension of the index 1 variables (must be > 0; if 0 then nSpecie_+1 is assumed)
        		0,                  //Dimension of the index 2 variables. Default nind2 = 0.
        		0,                  //Dimension of the index 2 variables. Default nind2 = 0.
        		1,                  //Switch for step size strategy: npred = 1 or 0  mod. predictive controller (Gustafsson) [default]; 
        		                    //                            If npred = 2  classical step size control (less safe; slightly faster)
        		0,                  //Default m1 = 0
        		0,                  //Default m2 = m1
        		false,              //If != 0 Jacobian matrix transformed to Hessenberg form 
        		                    //(advantageous for large systems with full Jacobian)
        		0,                  //Stopping criterion for Newton's method, usually chosen < 1. 
        		                    //Smaller values of fnewt make the code slower, but safer; Default min(0.03, sqrt(rtoler))
        		0.,                 //If quot1 < hnew/hold < quot2, then the step size is not changed.
        		0.,                 //This saves, together with a large thet, lu-decompositions and
                    				//computing time for large systems. for small systems one may have
                    				//quot1 = 1.0, quot2 = 1.2, for large full systems quot1 = 0.99,
                    				//quot2 = 2.0 might be good. Defaults quot1 = 1.0, quot2 = 1.2.

        		0.),                //Criterion for recomputing Jacobian; default 0.001; Increase thet to 0.1 if Jacobian costly
            // and solve ODEs for this cell

            Integrate();        

            // => solverX now contains result of ODE integration
            
            // do explicit underrelaxation of fine structures
            forAll(YStar_,specieI)
            {
                YStar_[specieI][cellI]=YStar_[specieI][cellI]*(1.-relaxFineStructures_) + solverX[specieI]*relaxFineStructures_;
            }
            
        } //end of cell loop

    } //end of direct integration (without table)

}

template<class CompType, class ThermoType>
Foam::scalar Foam::edcPSR<CompType, ThermoType>::solve
(
    const scalar startTime,
    const scalar deltaT
)
{   
    if(!this->chemistry_)
    {
        Info << "Chemistry disabled in const/chemistryProperties!!" << endl;
        return 0;
    }
    Info << "Selected chemistry model is EDC/PerfectlyStirredReactor" << endl;

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
    gammaStar_=pow(gammaL,2.);
    gammaStar_= min (0.7, gammaStar_);
    mDotStar_= 2.45*pow( epsilon/nu  , (1./2.) );


    // calculate fine structure composition according to fast chemistry approach
    Info << "Solving PSR equations..." << endl;
    updateYStar(startTime, deltaT);

    // if additional output requested, write TStar field

    if(additionalOutput_)
    {
        writeScalarField(gammaStar_,"gammaStar");
        writeScalarField(mDotStar_,"mDotStar");


        scalar TStarCell = 0.;
        scalar rhoStarCell = 0.;
        scalarField YStarCell(nSpecie_);
    
        scalarField TStar(rho.size());
        forAll(TStar,cellI)
        {
            dimensionedScalar hCell=this->thermo().hs().internalField()[cellI];
            scalar pCell=this->thermo().p().internalField()[cellI];
            scalar TCell=this->thermo().T().internalField()[cellI];
          
            forAll(Y_,specieI)
            {
                YStarCell[specieI]=YStar_[specieI][cellI];
            }
            iterateTAndRho(TStarCell, rhoStarCell, this->thermo().hs().internalField()[cellI],  pCell , YStarCell, TCell);
            TStar[cellI]=TStarCell;
        }
        writeScalarField(TStar,"TStar");
    }


    //assume X1*X2*X3=1 
    scalar X=1.;
    // calculate reaction rate coefficient                                          
    scalarField edcFactor = rho * gammaStar_ * mDotStar_ * X  / (1-gammaStar_*X);


                    
                    
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


//reimplement functions from RADAU5 library
template<class CompType, class ThermoType>
void Foam::edcPSR<CompType, ThermoType>::RADAU5derivative(double x, double *y, double *f) const
{
    // in case of PSR we do not use the value of x explicitly    
    
	Foam::scalarField vectorY(nSpecie_+1, 0.0);
	Foam::scalarField vectordYdt(nSpecie_+1, 0.0);
	
	// assign each value of array type "double" to a new array of type "scalar"
	forAll(vectorY,i)
	{
		vectorY[i] = y[i];
	}


    this->derivative(vectorY,vectordYdt);
	
	// and assign back from "scalar" to "double"
	forAll(vectorY,i)
	{
		f[i]=vectordYdt[i];
	}
	


}

// Jacobian left empty
template<class CompType, class ThermoType>
void Foam::edcPSR<CompType, ThermoType>::RADAU5jacobian(double x, double *y, double **J) const
{
    Info << endl << endl << endl << "************** JACOBIAN FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;
}

// Mass matrix is unity, therefore left empty
template<class CompType, class ThermoType>
void Foam::edcPSR<CompType, ThermoType>::RADAU5mass(double **M) const
{
    Info << endl << endl << endl << "************** MASS FUNCTION SHOULD NEVER BE CALLED ***************" << endl << endl<<endl;

} 


// ************************************************************************* //
