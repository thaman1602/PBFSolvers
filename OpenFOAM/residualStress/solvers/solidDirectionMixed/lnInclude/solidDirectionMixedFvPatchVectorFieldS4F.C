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

\*---------------------------------------------------------------------------*/

#include "solidDirectionMixedFvPatchVectorFieldS4F.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
//#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidDirectionMixedFvPatchVectorFieldS4F::solidDirectionMixedFvPatchVectorFieldS4F
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    limitCoeff_(1.0)
{
    // Lookup the solidModel object
    /*const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    if (solMod.solidModelDict().found("snGradLimitCoeff"))
    {
        limitCoeff_ =
            readScalar
            (
                solMod.solidModelDict().lookup("snGradLimitCoeff")
            );

        Info<< "snGradLimitCoeff: " << limitCoeff_ << endl;
    }*/
}


solidDirectionMixedFvPatchVectorFieldS4F::solidDirectionMixedFvPatchVectorFieldS4F
(
    const solidDirectionMixedFvPatchVectorFieldS4F& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    limitCoeff_(ptf.limitCoeff_)
{}


solidDirectionMixedFvPatchVectorFieldS4F::solidDirectionMixedFvPatchVectorFieldS4F
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF, dict),
    limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0))
{
    Info<< "Creating " << type() << " boundary condition" << endl;
    directionMixedFvPatchVectorField::evaluate();

    // Lookup the solidModel object
    /*const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    if (solMod.solidModelDict().found("snGradLimitCoeff"))
    {
        limitCoeff_ =
            readScalar
            (
                solMod.solidModelDict().lookup("snGradLimitCoeff")
            );

        Info<< "snGradLimitCoeff: " << limitCoeff_ << endl;
    }*/

    Info<< "Limiter coefficient: " << limitCoeff_ << endl;
}

solidDirectionMixedFvPatchVectorFieldS4F::solidDirectionMixedFvPatchVectorFieldS4F
(
    const solidDirectionMixedFvPatchVectorFieldS4F& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    limitCoeff_(ptf.limitCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidDirectionMixedFvPatchVectorFieldS4F::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidDirectionMixedFvPatchVectorFieldS4F::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void solidDirectionMixedFvPatchVectorFieldS4F::updateCoeffs()
{
    directionMixedFvPatchVectorField::updateCoeffs();
}

void solidDirectionMixedFvPatchVectorFieldS4F::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    bool secondOrder_(false);

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vcetors
    const vectorField k((I - sqr(n)) & delta);


    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Calc limited snGrad correction
    vectorField snGradCorrection
    (
      - (k & gradD.patchInternalField())*patch().deltaCoeffs()
    );

    if (limitCoeff_ < (1.0 - SMALL))
    {
        const vectorField uncorrectedSnGrad
        (
            (
               *this
             - patchInternalField()
            )*patch().deltaCoeffs()
        );

        const scalarField limiter
        (
            min
            (
                limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
               /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
                1.0
            )
        );

        snGradCorrection *= limiter;
    }

    const Field<vector> normalValue(transform(valueFraction(), refValue()));

    Field<vector> gradValue(patch().size());
    if (!secondOrder_)
    {
        gradValue =
            patchInternalField()
          - snGradCorrection/patch().deltaCoeffs()
          + refGrad()/patch().deltaCoeffs();
    }
    else
    {
        const vectorField nGradDP(n & gradD.patchInternalField());

        gradValue =
            patchInternalField()
          + (k & gradD.patchInternalField())
          + 0.5*(nGradDP + refGrad())/patch().deltaCoeffs();
    }

    const Field<vector> transformGradValue
    (
        transform(I - valueFraction(), gradValue)
    );

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}


Foam::tmp<Foam::Field<vector> >
solidDirectionMixedFvPatchVectorFieldS4F::snGrad() const
{
    const bool secondOrder_(false);

    const Field<vector> pif(this->patchInternalField());

    const Field<vector> normalValue
    (
        transform(this->valueFraction(), this->refValue())
    );

    // Unit normals
    const vectorField n(patch().nf());

    // Delta vectors
    const vectorField delta(patch().delta());

    // Correction vcetors
    const vectorField k((I - sqr(n)) & delta);

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    vectorField snGradCorrection
    (
      - (k & gradD.patchInternalField())
       *patch().deltaCoeffs()
    );

    if (limitCoeff_ < (1.0 - SMALL))
    {
        const vectorField uncorrectedSnGrad
        (
            (
               *this
              - patchInternalField()
            )*patch().deltaCoeffs()
        );

        const scalarField limiter
        (
            (
                min
                (
                    limitCoeff_*mag(uncorrectedSnGrad + snGradCorrection)
                   /((1 - limitCoeff_)*mag(snGradCorrection) + SMALL),
                    1.0
                )
            )
        );

        snGradCorrection *= limiter;
    }

    Field<vector> gradValue
    (
        pif
      - snGradCorrection/patch().deltaCoeffs()
      + refGrad()/patch().deltaCoeffs()
    );

    if (secondOrder_)
    {
        const vectorField nGradDP(n & gradD.patchInternalField());

        gradValue =
            patchInternalField()
          + (k & gradD.patchInternalField())
          + 0.5*(nGradDP + refGrad())/patch().deltaCoeffs();
    }

    Field<vector> transformGradValue
    (
        transform(I - this->valueFraction(), gradValue)
    );

    if (secondOrder_)
    {
        const vectorField nGradDP(n & gradD.patchInternalField());

        return
            2.0
           *(
               normalValue + transformGradValue
             - (pif + (k&gradD.patchInternalField()))
           )*patch().deltaCoeffs()
         - nGradDP;
    }

    return
    (
        normalValue + transformGradValue
      - (
            pif
          - snGradCorrection/this->patch().deltaCoeffs()
        )
    )*patch().deltaCoeffs();
}

// Write
void solidDirectionMixedFvPatchVectorFieldS4F::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    os.writeKeyword("limitCoeff")
        << limitCoeff_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidDirectionMixedFvPatchVectorFieldS4F);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
