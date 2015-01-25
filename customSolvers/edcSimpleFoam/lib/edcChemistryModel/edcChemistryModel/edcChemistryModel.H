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
    Foam::edcChemistryModel

Description
    Class to calculate chemical reaction rates with Eddy-Dissipation Concept (EDC)
    For details about EDC theory please refer to:
        B. Magnussen, "The Eddy-Dissipation Concept: 
        A bridge between science and technology" (2005)

    OpenFoam implementation by
        Bjorn Lilleberg and Dominik Christ

\*---------------------------------------------------------------------------*/

#ifndef edcChemistryModel_H
#define edcChemistryModel_H

#include "basicChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hsReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class edcChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class edcChemistryModel
:
    public basicChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        edcChemistryModel(const edcChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const edcChemistryModel&);


protected:

    // Protected data

        //- Thermo package
//        autoPtr<hCombustionThermo> thermo_;
        autoPtr<hsReactionThermo> thermo_;


public:

    //- Runtime type information
    TypeName("edcChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        edcChemistryModel,
        fvMesh,
        (
            const fvMesh& mesh,
            const word& compTypeName,
            const word& thermoTypeName
        ),
        (mesh, compTypeName, thermoTypeName)
    );


    // Constructors

        //- Construct from mesh
        edcChemistryModel(const fvMesh& mesh, const word& thermoTypeName);


    //- Selector
    static autoPtr<edcChemistryModel> New(const fvMesh& mesh);


    //- Destructor
    virtual ~edcChemistryModel();

    // Member Functions

        //- Return access to the thermo package
//        inline hCombustionThermo& thermo();
        inline hsReactionThermo& thermo();

        //- Return const access to the thermo package
//        inline const hCombustionThermo& thermo() const;
        inline const hsReactionThermo& thermo() const;        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edcChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
