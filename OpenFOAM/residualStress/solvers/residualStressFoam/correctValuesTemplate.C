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

\*---------------------------------------------------------------------------*/
#include "emptyFvPatchFields.H"
#include "correctValues.H"

namespace Foam
{

template<class Type>
void correctValues
(
   GeometricField<Type, fvPatchField, volMesh>& field,
   const word& interfaceName
)
{
    // interface patch ID in local mesh
    label localSolidInterfaceID = field.mesh().boundaryMesh().findPatchID(interfaceName);   
    // fvPatchSclarField associoated with interface
    fvPatchField<Type>& pField = field.boundaryField()[localSolidInterfaceID];

    const Field<Type> pInternalField = pField.patchInternalField();

    if (pField.size() != 0)
    {
        forAll (pInternalField, cellI)
        {
            pField[cellI] = pInternalField[cellI];
        }
    }
}

template<class Type>
void correctValues
(
   GeometricField<Type, fvsPatchField, surfaceMesh>& field,
   const word& interfaceName
)
{
    // interface patch ID in local mesh
    label localSolidInterfaceID = field.mesh().boundaryMesh().findPatchID(interfaceName);   
    // fvPatchSclarField associoated with interface
    fvsPatchField<Type>& pField = field.boundaryField()[localSolidInterfaceID];

    forAll (pField, cellI)
    {
      //Info << "Not corrected[" << cellI << "]: patch: " << pField[cellI] << endl; 
      //pField[cellI] = pInternalField[cellI];
      //Info << "Corrected[" << cellI << "]:  patch: " << pField[cellI] << endl; 
    }
}

} //end of namespace Foam

