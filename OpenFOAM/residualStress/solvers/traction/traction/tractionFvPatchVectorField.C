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

\*---------------------------------------------------------------------------*/

#include "tractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tractionFvPatchVectorField::
tractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(false),
    setEffectiveTraction_(false),
    limitCoeff_(1.0),
    relaxFac_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


tractionFvPatchVectorField::
tractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    setEffectiveTraction_
    (
        dict.lookupOrDefault<Switch>("setEffectiveTraction", false)
    ),
    limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0)),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check if traction is time-varying
    if (dict.found("tractionSeries"))
    {
        Info<< "    traction is time-varying" << endl;
        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        Info<< "    pressure is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

    if (secondOrder_)
    {
        Info<< "    second order correction" << endl;
    }

    if (setEffectiveTraction_)
    {
        Info<< "    set effective traction" << endl;
    }

    if (limitCoeff_)
    {
        Info<< "    limiter coefficient: " << limitCoeff_ << endl;
    }

    if (relaxFac_ < 1.0)
    {
        Info<< "    relaxation factor: " << relaxFac_ << endl;
    }
}


tractionFvPatchVectorField::
tractionFvPatchVectorField
(
    const tractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    traction_(mapper(stpvf.traction_)),
    pressure_(mapper(stpvf.pressure_)),
#else
    traction_(stpvf.traction_, mapper),
    pressure_(stpvf.pressure_, mapper),
#endif
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    secondOrder_(stpvf.secondOrder_),
    setEffectiveTraction_(stpvf.setEffectiveTraction_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


tractionFvPatchVectorField::
tractionFvPatchVectorField
(
    const tractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    secondOrder_(stpvf.secondOrder_),
    setEffectiveTraction_(stpvf.setEffectiveTraction_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


tractionFvPatchVectorField::
tractionFvPatchVectorField
(
    const tractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    secondOrder_(stpvf.secondOrder_),
    setEffectiveTraction_(stpvf.setEffectiveTraction_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void tractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(traction_, traction_);
    m(pressure_, pressure_);
#else
    traction_.autoMap(m);
    pressure_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void tractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const tractionFvPatchVectorField& dmptf =
        refCast<const tractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void tractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    scalarField press = pressure_;
    if (setEffectiveTraction_)
    {
        const fvPatchField<scalar>& p =
            patch().lookupPatchField<volScalarField, scalar>("p");

        // Remove the dynamic pressure component: this will force the effective
        // traction to be enforced rather than the total traction
        press -= p;
    }

    // solids4foam Sep-21
    // // Lookup the solidModel object
    // const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    // // Set surface-normal gradient on the patch corresponding to the desired
    // // traction
    // gradient() =
    //     relaxFac_*solMod.tractionBoundarySnGrad
    //     (
    //         traction_, press, patch()
    //     )
    //   + (1.0 - relaxFac_)*gradient();

    // Method which direclty looks sigma and impK fields and assumes small
    // strains

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

    // Return patch snGrad
    const vectorField newGradient =
    (
        (traction_ - n*pressure_)
      - (n & (pSigma - impK*pGradD))
    )/impK;

    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    gradient() = relaxFac_*newGradient + (1.0 - relaxFac_)*gradient();

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void tractionFvPatchVectorField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
#ifdef OPENFOAMESIORFOUNDATION
            "grad(" + internalField().name() + ")"
#else
            "grad(" + dimensionedInternalField().name() + ")"
#endif
        );

    // Face unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Non-orthogonal correction vectors
    const vectorField k = ((I - sqr(n)) & delta);

    if (secondOrder_)
    {
        const vectorField dUP = (k & gradField.patchInternalField());
        const vectorField nGradUP = (n & gradField.patchInternalField());

        Field<vector>::operator=
        (
            patchInternalField()
          + dUP
          + 0.5*(gradient() + nGradUP)/patch().deltaCoeffs()
        );
    }
    else
    {

        Field<vector>::operator=
        (
            patchInternalField()
          + (k & gradField.patchInternalField())
          + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


void tractionFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //fixedGradientFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "traction", traction_);
#else
        traction_.writeEntry("traction", os);
#endif
    }

    if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "pressure", pressure_);
#else
        pressure_.writeEntry("pressure", os);
#endif
    }
    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
    os.writeKeyword("setEffectiveTraction")
        << setEffectiveTraction_ << token::END_STATEMENT << nl;
    os.writeKeyword("limitCoeff")
        << limitCoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << nl;

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
    writeEntry(os, "gradient", gradient());
#else
    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, tractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
