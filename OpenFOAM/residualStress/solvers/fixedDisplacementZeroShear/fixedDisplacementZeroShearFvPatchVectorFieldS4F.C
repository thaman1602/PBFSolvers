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

#include "fixedDisplacementZeroShearFvPatchVectorFieldS4F.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorFieldS4F::
fixedDisplacementZeroShearFvPatchVectorFieldS4F
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorFieldS4F(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_(),
    forceZeroShearGrad_(false),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0)
{}


fixedDisplacementZeroShearFvPatchVectorFieldS4F::
fixedDisplacementZeroShearFvPatchVectorFieldS4F
(
    const fixedDisplacementZeroShearFvPatchVectorFieldS4F& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorFieldS4F(ptf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    totalDisp_(mapper(ptf.totalDisp_)),
#else
    totalDisp_(ptf.totalDisp_, mapper),
#endif
    dispSeries_(ptf.dispSeries_),
    forceZeroShearGrad_(ptf.forceZeroShearGrad_),
    traction_(ptf.traction_, mapper),
    pressure_(ptf.pressure_, mapper)
{}


fixedDisplacementZeroShearFvPatchVectorFieldS4F::
fixedDisplacementZeroShearFvPatchVectorFieldS4F
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorFieldS4F(p, iF),
    totalDisp_("value", dict, p.size()),
    dispSeries_(),
    forceZeroShearGrad_
    (
        dict.lookupOrDefault<Switch>("forceZeroShearGrad", false)
    ),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size())
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        Info<< "    displacement is time-varying" << endl;
        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        refValue() = dispSeries_(this->db().time().timeOutputValue());
    }
    else if (dict.found("value"))
    {
        refValue() = vectorField("value", dict, p.size());
    }
    else
    {
        FatalErrorIn
        (
            "fixedDisplacementZeroShearFvPatchVectorFieldS4F::"
            "fixedDisplacementZeroShearFvPatchVectorFieldS4F"
        )   << "value entry not found for patch " << patch().name()
            << abort(FatalError);
    }

    this->refGrad() = vector::zero;

    this->valueFraction() = sqr(patch().nf());

    Field<vector> normalValue
    (
        transform(valueFraction(), refValue())
    );

    Field<vector> gradValue
    (
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs()
    );

    Field<vector> transformGradValue
    (
        transform(I - valueFraction(), gradValue)
    );

    Field<vector>::operator=(normalValue + transformGradValue);
}


fixedDisplacementZeroShearFvPatchVectorFieldS4F::
fixedDisplacementZeroShearFvPatchVectorFieldS4F
(
    const fixedDisplacementZeroShearFvPatchVectorFieldS4F& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorFieldS4F(ptf, iF),
    totalDisp_(ptf.totalDisp_),
    dispSeries_(ptf.dispSeries_),
    forceZeroShearGrad_(ptf.forceZeroShearGrad_),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorFieldS4F::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorFieldS4F::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(totalDisp_, totalDisp_);;
#else
    totalDisp_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorFieldS4F::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorFieldS4F::rmap(ptf, addr);

    const fixedDisplacementZeroShearFvPatchVectorFieldS4F& dmptf =
        refCast<const fixedDisplacementZeroShearFvPatchVectorFieldS4F>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedDisplacementZeroShearFvPatchVectorFieldS4F::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    vectorField disp = totalDisp_;

    if (dispSeries_.size())
    {
        disp = dispSeries_(this->db().time().timeOutputValue());
    }

#ifdef OPENFOAMESIORFOUNDATION
    if (internalField().name() == "D")
#else
    if (dimensionedInternalField().name() == "D")
#endif
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    // Set displacement
    refValue() = disp;

    // Set gradient to zero to force zero shear traction
    if (forceZeroShearGrad_)
    {
        refGrad() = vector::zero;
    }
    else
    {
        // Calculate the shear gradient such that the shear traction is zero

        // Patch mechanical property
        const scalarField& impK =
            patch().lookupPatchField<volScalarField, scalar>("impK");

        // Patch gradient of displacement
        const tensorField& pGradD =
            patch().lookupPatchField<volTensorField, tensor>("grad(D)");

        // Patch stress
        const symmTensorField& pSigma =
            patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");

        // Patch unit normals
        const vectorField n = patch().nf();

        // Set gradient to force zero shear traction
        refGrad() =
                    (
                        (traction_ - n*pressure_)
                        - (n & (pSigma - impK*pGradD))
                    )/impK;
    }

    solidDirectionMixedFvPatchVectorFieldS4F::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearFvPatchVectorFieldS4F::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    os.writeKeyword("forceZeroShearGrad")
        << forceZeroShearGrad_ << token::END_STATEMENT << nl;

    solidDirectionMixedFvPatchVectorFieldS4F::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearFvPatchVectorFieldS4F
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
