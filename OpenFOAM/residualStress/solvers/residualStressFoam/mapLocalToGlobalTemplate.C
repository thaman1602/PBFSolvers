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
#include "fvCFD.H"
#include "emptyFvPatchFields.H"
#include "correctValues.H"

namespace Foam
{

  /********** MAP VALUES **********/
template<class Type>
void mapLocalToGlobal
(
 GeometricField<Type, fvPatchField, volMesh>& global,
 const GeometricField<Type, fvPatchField, volMesh>& local,
 const fvMeshSubset& subsetMesh
)
{
  const fvMesh& localMesh = local.mesh();
  const fvMesh& globalMesh = global.mesh();
//at the end of time step
//return temperature field from local to global
//list relating cell centre from local to global mehs
const labelList& cellMap = subsetMesh.cellMap();
forAll(local, cellI)
{
  //return ID for cell of global mesh
  const label baseMeshCellID = cellMap[cellI];
  global[baseMeshCellID] = local[cellI];
}

/**** Mapping face values ****/
//connect local and global patch ID's
const labelList& patchMap = subsetMesh.patchMap();
//connect local and global face ID's - internal faces included
const labelList& faceMap = subsetMesh.faceMap();

forAll (localMesh.boundary(), patchI)
{
  //current patch
  const fvPatch& curPatch = localMesh.boundary()[patchI];
  //attached global-patch ID
  const label baseMeshPatchID = patchMap[patchI];

  if (curPatch.type() != "empty")
  {
    //starting face ID for current patch in global
    // addresing
    const label start = curPatch.patch().start();

    //temperature sub-Field
    const Field<Type>& pLocal = local.boundaryField()[patchI];

    forAll (pLocal, faceI)
    {
      // get appropiate global mesh face ID (global addresing)
      // with start all internal fields and previous patches
      // are skiped
      label baseMeshFaceGlobalID = faceMap[start + faceI];
      if (baseMeshFaceGlobalID >= globalMesh.nInternalFaces())
      {
        label globalStart = globalMesh.boundary()[baseMeshPatchID].patch().start();
        label baseMeshFaceID = baseMeshFaceGlobalID - globalStart;

        global.boundaryField()[baseMeshPatchID][baseMeshFaceID] = local.boundaryField()[patchI][faceI];
      }
    }
  }
}
}

} //end of namespace Foam

