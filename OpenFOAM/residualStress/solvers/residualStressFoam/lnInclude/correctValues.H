/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    correctValues()
    Patch on the interface does not contain faces in the base mesh and 
    in submesh random values are put on those faces. This function puts
    adjoining cell values on the interface patch faces. correctInterface()
    is used just for the field with specified BC's, correctValues() has to be
    used for every other field.

    checkInterfacePatchName(mesh)
    checks is patch under "solidInterface" is found - for now that name is
    hardcoded in the curing library

SourceFiles

\*---------------------------------------------------------------------------*/
#ifndef correctInterface_H
#define correctInterface_H

#include "fvCFD.H"
#include "FvMeshSubset.H"

namespace Foam
{

void checkInterfacePatchName
(
   const fvMesh& mesh
);

template<class Type>
void correctValues
(
    GeometricField<Type, fvPatchField, volMesh>& field,
   const word& interfaceName
);

template<class Type>
void correctValues
(
   GeometricField<Type, fvsPatchField, surfaceMesh>& field,
   const word& interfaceName
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "correctValuesTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
