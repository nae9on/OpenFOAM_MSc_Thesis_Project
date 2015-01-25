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

InClass
    Foam::edcChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"
#include "thermoPhysicsTypes.H"

#include "edcChemistryModel.H"
#include "edm.H"
#include "edmKinetic.H"
#include "edcFastChem.H"
#include "edcLE.H"
#include "edcPSR.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeChemistryModel
    (
        edm,
        edcChemistryModel,
        gasThermoPhysics
    );

    makeChemistryModel
    (
        edmKinetic,
        edcChemistryModel,
        gasThermoPhysics
    );    

    makeChemistryModel
    (
        edcFastChem,
        edcChemistryModel,
        gasThermoPhysics
    ); 
    
    makeChemistryModel
    (
        edcLE,
        edcChemistryModel,
        gasThermoPhysics
    );     
    makeChemistryModel
    (
        edcPSR,
        edcChemistryModel,
        gasThermoPhysics
    );     

    
}

// ************************************************************************* //
