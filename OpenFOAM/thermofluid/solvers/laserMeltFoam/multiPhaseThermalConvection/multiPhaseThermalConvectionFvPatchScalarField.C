/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "multiPhaseThermalConvectionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiPhaseThermalConvectionFvPatchScalarField::multiPhaseThermalConvectionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    KName_("undefined"),
    hc_(p.size(), 0),
    Tinf_("Tinf", dimTemperature, 293)
{}


Foam::multiPhaseThermalConvectionFvPatchScalarField::multiPhaseThermalConvectionFvPatchScalarField
(
    const multiPhaseThermalConvectionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    KName_(ptf.KName_),
    hc_(ptf.hc_, mapper),
    Tinf_(ptf.Tinf_)
{}


Foam::multiPhaseThermalConvectionFvPatchScalarField::multiPhaseThermalConvectionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    KName_(dict.lookupOrDefault<word>("thermalConductivityName", "k")),
    hc_("hc", dict, p.size()),
    Tinf_("Tinf", dimTemperature, readScalar(dict.lookup("Tinf")))
{
    Info<< patch().name() << ": multiPhaseThermalConvection" << endl;

    fvPatchField<scalar>::operator=(patchInternalField());
}


Foam::multiPhaseThermalConvectionFvPatchScalarField::multiPhaseThermalConvectionFvPatchScalarField
(
    const multiPhaseThermalConvectionFvPatchScalarField& wbppsf
)
:
    fixedValueFvPatchScalarField(wbppsf),
    KName_(wbppsf.KName_),
    hc_(wbppsf.hc_),
    Tinf_(wbppsf.Tinf_)
{}


Foam::multiPhaseThermalConvectionFvPatchScalarField::multiPhaseThermalConvectionFvPatchScalarField
(
    const multiPhaseThermalConvectionFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wbppsf, iF),
    KName_(wbppsf.KName_),
    hc_(wbppsf.hc_),
    Tinf_(wbppsf.Tinf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::multiPhaseThermalConvectionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::multiPhaseThermalConvectionFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const multiPhaseThermalConvectionFvPatchScalarField& tiptf =
        refCast<const multiPhaseThermalConvectionFvPatchScalarField>(ptf);

    //alpha_.resize(tiptf.alpha_.size());
    hc_.rmap(tiptf.hc_, addr);
}


void Foam::multiPhaseThermalConvectionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalarField delta = 1.0/this->patch().deltaCoeffs() + SMALL;
    scalarField TP = this->patchInternalField();

if (this->db().foundObject<surfaceScalarField>(KName_))
{
    // Lookup thermal diffusivity i.e. conductivity
    const fvsPatchField<scalar>& K =
        patch().lookupPatchField<surfaceScalarField, scalar>(KName_);

    Field<scalar>::operator=
    (
        (hc_*Tinf_.value() + K*TP/delta)
       /(K/delta + hc_ + SMALL)
    );
    fvPatchField<scalar>::evaluate();
} else
{
    // Lookup thermal diffusivity i.e. conductivity
         const fvPatchField<scalar>& K =
                 patch().lookupPatchField<volScalarField, scalar>(KName_);
    
                     Field<scalar>::operator=
                         (
                                 (hc_*Tinf_.value() + K*TP/delta)
                                        /(K/delta + hc_ + SMALL)
                            );
 }
   fvPatchField<scalar>::evaluate();
}


/*Foam::tmp<Foam::Field<Foam::scalar> >
Foam::multiPhaseThermalConvectionFvPatchScalarField::snGrad() const
{
    scalarField delta = 1.0/this->patch().deltaCoeffs() + SMALL;
    scalarField TP = this->patchInternalField();

    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

//     scalarField gradient =
//          DT*alpha_*(Tinf_.value() - TP)
//         /(DT + alpha_*delta + SMALL);

//     scalarField gradient =
//          alpha_*(Tinf_.value() - TP)
//         /(DT + alpha_*delta + SMALL);

    return tmp<Field<scalar> >
    (
         alpha_*(Tinf_.value() - TP)
        /(DT + alpha_*delta + SMALL)
    );
}*/


/*Foam::tmp<Foam::Field<Foam::scalar> >
Foam::multiPhaseThermalConvectionFvPatchScalarField::gradientInternalCoeffs() const
{
    scalarField delta = 1.0/this->patch().deltaCoeffs() + SMALL;

    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

//     return tmp<Field<scalar> >
//     (
//         -pTraits<scalar>::one*DT*alpha_/(DT + alpha_*delta + SMALL)
//     );

    return tmp<Field<scalar> >
    (
        -pTraits<scalar>::one*alpha_/(DT + alpha_*delta + SMALL)
    );
}


Foam::tmp<Foam::Field<Foam::scalar> >
Foam::multiPhaseThermalConvectionFvPatchScalarField::gradientBoundaryCoeffs() const
{
    scalarField delta = 1.0/this->patch().deltaCoeffs() + SMALL;

    const fvPatchField<scalar>& DT =
        patch().lookupPatchField<volScalarField, scalar>(DTName_);

    return alpha_*Tinf_.value()/(DT + alpha_*delta + SMALL);
}*/


void Foam::multiPhaseThermalConvectionFvPatchScalarField::write(Ostream& os) const
{
    //fvPatchScalarField::write(os);
    //writeEntry("value", os);

    os.writeKeyword("thermalConductivityName")
        << KName_ << token::END_STATEMENT << nl;
    hc_.writeEntry("hc", os);
    os.writeKeyword("Tinf") << Tinf_.value() << token::END_STATEMENT << nl;

    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        multiPhaseThermalConvectionFvPatchScalarField
    );
}

// ************************************************************************* //
