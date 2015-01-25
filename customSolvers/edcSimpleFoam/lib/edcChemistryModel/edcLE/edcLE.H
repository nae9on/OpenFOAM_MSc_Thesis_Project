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
    Foam::edcLE

Description
    Implements edcChemistryModel with reaction rates from Eddy-Dissipation Concept 
    using Fast Chemistry approach and Local Extinction database

SourceFiles
    edcLEI.H
    edcLE.C

\*---------------------------------------------------------------------------*/

#ifndef edcLE_H
#define edcLE_H

#include "hsReactionThermo.H"
#include "Reaction.H"
#include "interpolateXY.H"
#include "itoa.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class edcLE Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class edcLE
:
    public CompType

{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const edcLE&);


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
        
        // under-relaxation factor for fine structure mass fractions
        scalar relaxFineStructures_;

        // name of oxidiser species for each reaction
        label indexOxidiser_;

        // name of main fuel species used for calculating Local Extinction
        label indexMainFuel_;

        // name of main fuel species used for calculating Local Extinction
        label indexMainReaction_;        
        
        // stoichiometric fuel/oxidiser ratio
        scalar fStoich_;
        
        // Local Extinction coefficients
        // temperature level of given curve
        scalarField LETemp_;
        
        // smallest chemical time scale
        scalar autoIgnitionTemperature_;

        // smallest chemical time scale
        scalarField tauChMin_;
        
       
        // is additional output required?
        bool additionalOutput_;
        
        scalarField tauCh_;

        scalar reason_;

    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<scalarField>& RR();
        

        PtrList<scalarField> calculateYStar();

        void writeScalarField(scalarField X, word fieldName);
        
        bool isLocalExtinction(scalar temp, scalar equivalenceRatio, scalar tauEdc);

public:

    //- Runtime type information
    TypeName("edcLE");


    // Constructors

        //- Construct from components
        edcLE
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        );


    //- Destructor
    virtual ~edcLE();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        // Chemistry model functions (implementing abstract functions in
        // chemistryModel.H)

            //- Return const access to the chemical source terms
            inline tmp<volScalarField> RR(const label i) const;

            //- Solve the reaction system for the given start time and time
            //  step and return the characteristic time
            virtual scalar solve(const scalar t0, const scalar deltaT);

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

#include "edcLEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "edcLE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
