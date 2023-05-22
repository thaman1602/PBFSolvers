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
    fixedDisplacementZeroShearFvPatchVectorFieldS4F

Description
    The component of the displacement normal to the patch is applied,
    and the shear traction is set to zero.

    Patch contains non-orthogonal correction.

    Optionally time-varying by specifying a time series.

SourceFiles
    fixedDisplacementZeroShearFvPatchVectorFieldS4F.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.
    Based on velocityDisplacement by Aleksandar Karac

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementZeroShearFvPatchVectorFieldS4F_H
#define fixedDisplacementZeroShearFvPatchVectorFieldS4F_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorFieldS4F.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementZeroShearFvPatchVectorFieldS4F
:
    public solidDirectionMixedFvPatchVectorFieldS4F
{
    // Private data

        //- Constant total displacement field
        vectorField totalDisp_;

        //- Displacement time series
        interpolationTable<vector> dispSeries_;

        //- Switch to force the shear gradient to be zero
        const Switch forceZeroShearGrad_;

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;


public:

    //- Runtime type information
    TypeName("fixedDisplacementZeroShearS4F");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementZeroShearFvPatchVectorFieldS4F
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementZeroShearFvPatchVectorFieldS4F
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedDisplacementZeroShearFvPatchVectorFieldS4F
        //  onto a new patch
        fixedDisplacementZeroShearFvPatchVectorFieldS4F
        (
            const fixedDisplacementZeroShearFvPatchVectorFieldS4F&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedDisplacementZeroShearFvPatchVectorFieldS4F(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedDisplacementZeroShearFvPatchVectorFieldS4F
        (
            const fixedDisplacementZeroShearFvPatchVectorFieldS4F&,
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
                new fixedDisplacementZeroShearFvPatchVectorFieldS4F(*this, iF)
            );
        }

    // Destructor

        virtual ~fixedDisplacementZeroShearFvPatchVectorFieldS4F()
        {}


    // Member functions

            virtual const vectorField& traction() const
            {
                return traction_;
            }

            virtual vectorField& traction()
            {
                return traction_;
            }

            virtual const scalarField& pressure() const
            {
                return pressure_;
            }

            virtual scalarField& pressure()
            {
                return pressure_;
            }

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


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
