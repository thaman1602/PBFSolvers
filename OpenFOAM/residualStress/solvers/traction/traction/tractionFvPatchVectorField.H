/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    tractionFvPatchVectorField

Description
    Fixed traction boundary condition for solid solvers.

    This is a modified version of solidTraction from solids4foam such that it no
    longer requires other classes within solids4foam.

Author
    Philip Cardiff, UCD.
    Zeljko Tukovic, FSB Zagreb.

SourceFiles
    solidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionFvPatchVectorField_H
#define tractionFvPatchVectorField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class tractionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        //- Traction
        vectorField traction_;

        //- Pressure
        scalarField pressure_;

        //- Traction time series
        interpolationTable<vector> tractionSeries_;

        //- Pressure time series
        interpolationTable<scalar> pressureSeries_;

        //- Second order accuracy correction
        const Switch secondOrder_;

        //- Switch so that the effective traction is enforced rather than the
        //  total traction. This only applies to porous flow based solids
        const Switch setEffectiveTraction_;

        //- Limiter coefficient for snGrad correction
        const scalar limitCoeff_;

        //- Optional: under-relaxation factor
        const scalar relaxFac_;


public:

    //- Runtime type information
    TypeName("traction");


    // Constructors

        //- Construct from patch and internal field
        tractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  tractionFvPatchVectorField onto a new patch
        tractionFvPatchVectorField
        (
            const tractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionFvPatchVectorField
        (
            const tractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionFvPatchVectorField
        (
            const tractionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access

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
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
#ifdef OPENFOAMESIORFOUNDATION
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
#else
            const Pstream::commsTypes commsType = Pstream::blocking
#endif
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
