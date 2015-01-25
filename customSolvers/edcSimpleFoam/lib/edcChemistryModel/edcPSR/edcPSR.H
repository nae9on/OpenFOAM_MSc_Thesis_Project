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

Class
    Foam::edcPSR

Description
    Implements edcChemistryModel with reaction rates from Eddy-Dissipation Concept 
    using Fast Chemistry approach for fine structure mass fractions

SourceFiles
    edcPSRI.H
    edcPSR.C

\*---------------------------------------------------------------------------*/

#ifndef edcPSR_H
#define edcPSR_H

#include "hsReactionThermo.H"
#include "Reaction.H"
#include "binaryTree.H"
#include "StiffIntegratorT.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class edcPSR Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class edcPSR
:
    public CompType,
    public StiffIntegratorT

{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const edcPSR&);


protected:

    // Private data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction<ThermoType> >& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;
        
        //- Chemical source term
        PtrList<scalarField> RR_;
        
        // fine structure mass fractions
        PtrList<scalarField> YStar_;

        // mass transfer rate to and from fine structures
        scalarField mDotStar_;        

        // mass fraction of fine structures in control volume
        scalarField gammaStar_;    
        
        // under-relaxation factor for fine structure mass fractions
        scalar relaxFineStructures_;

        // is additional output required?
        bool additionalOutput_;
        
        // use binary tree?
        bool useBinaryTree_;       
        
        // binary tree to store results of YStar integrations
        binaryTree* resultTable_;
        
        // max number of table entries
        label tableSize_;

        // allowable error for binary tree lookups
        scalar tableErr_;        
        
        //normalization values for resultTable
        scalar hsMax_;
        scalar gammaStarMax_;
        scalar mDotStarMax_;
        
        // parameters for RADAU5 solver
        scalar absoluteTolerance_;
        scalar relativeTolerance_;
        scalar maxIterations_;  
        
        // pointer to RADAU5 solution vector
        double* solverX;        
    
        //temperature limits
        scalar TMin_;
        scalar TMax_;
        


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<scalarField>& RR();
        
        void updateYStar(scalar startTime, scalar stopTime);

        void writeScalarField(scalarField X, word fieldName);

        void iterateTAndRho(scalar& TResult, scalar& rho, scalar h, scalar p, scalarField& Y, scalar TInit) const;

        void derivative(scalarField&  YStarCell, scalarField& dYdtCell) const;
        
public:

    //- Runtime type information
    TypeName("edcPSR");


    // Constructors

        //- Construct from components
        edcPSR
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        );


    //- Destructor
    virtual ~edcPSR();



    // Member Functions

    // implement functions from RADAU5 solver
        void RADAU5derivative(double x, double *y, double *f) const;
        void RADAU5jacobian(double x, double *y, double **J) const;
        void RADAU5mass(double **M) const;

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        //- dc/dt = omega, rate of change in concentration, for each species
        scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- return the reaction rate for reaction r and the reference
        //  species and charateristic times
        scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        // Chemistry model functions (implementing abstract functions in
        // chemistryModel.H)

            //- Return const access to the chemical source terms
            inline tmp<volScalarField> RR(const label i) const;

            //- Solve the reaction system for the given start time and time
            //  step and return the characteristic time
            virtual scalar solve(const scalar startTime, const scalar deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release
            virtual tmp<volScalarField> dQ() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edcPSRI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "edcPSR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
