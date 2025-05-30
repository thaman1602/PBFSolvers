/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidDirectionMixedFvPatchVectorFieldS4F

Description
    directionMixed with non-orthogonal correction for the diffusion term

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

SourceFiles
    solidDirectionMixedFvPatchVectorFieldS4F.C

\*---------------------------------------------------------------------------*/

#ifndef solidDirectionMixedFvPatchVectorFieldS4F_H
#define solidDirectionMixedFvPatchVectorFieldS4F_H

#include "directionMixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidDirectionMixedFvPatchVectorFieldS4F
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Limiter coefficient for snGrad correction
        scalar limitCoeff_;


public:

    //- Runtime type information
    TypeName("solidDirectionMixed");


    // Constructors

        //- Construct from patch and internal field
        solidDirectionMixedFvPatchVectorFieldS4F
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidDirectionMixedFvPatchVectorFieldS4F
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidDirectionMixedFvPatchVectorFieldS4F
        //  onto a new patch
        solidDirectionMixedFvPatchVectorFieldS4F
        (
            const solidDirectionMixedFvPatchVectorFieldS4F&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidDirectionMixedFvPatchVectorFieldS4F(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidDirectionMixedFvPatchVectorFieldS4F
        (
            const solidDirectionMixedFvPatchVectorFieldS4F&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidDirectionMixedFvPatchVectorFieldS4F(*this, iF)
            );
        }

    // Destructor

        virtual ~solidDirectionMixedFvPatchVectorFieldS4F()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- evaluate the patch
            virtual void evaluate(const Pstream::commsTypes);

            //- snGrad with non-orthogonal correction
            tmp<Foam::Field<vector> > snGrad() const;


        //- Write

            virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
