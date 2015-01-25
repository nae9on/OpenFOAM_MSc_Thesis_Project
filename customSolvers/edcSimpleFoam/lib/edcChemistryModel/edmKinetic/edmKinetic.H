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
    Foam::edmKinetic

Description
    Implements edcChemistryModel with reaction rates from Eddy-Dissipation Model (i.e. fast chemistry) approach

SourceFiles
    edmKineticI.H
    edmKinetic.C

\*---------------------------------------------------------------------------*/

#ifndef edmKinetic_H
#define edmKinetic_H

#include "hsReactionThermo.H"
#include "Reaction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class edmKinetic Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class edmKinetic
:
    public CompType

{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const edmKinetic&);


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

        //- Chemistry solver
//        autoPtr<chemistrySolver<CompType, ThermoType> > solver_;

        //- Chemical source term
        PtrList<scalarField> RR_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<scalarField>& RR();


public:

    //- Runtime type information
    TypeName("edmKinetic");


    // Constructors

        //- Construct from components
        edmKinetic
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        );


    //- Destructor
    virtual ~edmKinetic();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
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

#include "edmKineticI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "edmKinetic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
