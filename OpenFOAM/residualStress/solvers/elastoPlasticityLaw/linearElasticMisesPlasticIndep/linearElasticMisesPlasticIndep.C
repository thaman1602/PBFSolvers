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

#include "linearElasticMisesPlasticIndep.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMisesPlasticIndep, 0);

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar linearElasticMisesPlasticIndep::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label linearElasticMisesPlasticIndep::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar linearElasticMisesPlasticIndep::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar linearElasticMisesPlasticIndep::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::linearElasticMisesPlasticIndep::updatePlasticity
(
    symmTensor& plasticN,          // Plastic return direction
    scalar& DLambda,               // Plastic multiplier increment
    scalar& DSigmaY,               // Increment of yield stress
    scalar& sigmaY,                // Yield stress
    const scalar sigmaYOld,        // Yield stress old time
    const scalar fTrial,           // Trial yield function
    const symmTensor& sTrial,      // Trial deviatoric stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagBE          // Max strain increment magnitude
) const
{
    // Calculate DLambda/DEpsilonPEq
    if (fTrial < SMALL)
    {
        // Elasticity
        plasticN = symmTensor(I);
        DLambda = 0.0;
        DSigmaY = 0.0;
        sigmaY = sigmaYOld;
    }
    else
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrial);
        if (magS > SMALL)
        {
            plasticN = sTrial/magS;
        }
        else
        {
            // Deviatoric stress is zero so plasticN value does not matter, but
            // we will set it to the identity
            plasticN = symmTensor(I);
        }

        if (nonLinearPlasticity_)
        {
            // Update plastic multiplier (DLambda) and current yield stress
            // (sigmaY)
            newtonLoop
            (
                DLambda,
                sigmaY,
                epsilonPEqOld,
                magS,
                mu_.value(),
                maxMagBE
            );

            // Update increment of yield stress
            DSigmaY = sigmaY - sigmaYOld;
        }
        else
        {
            // Update DLambda
            DLambda = fTrial/(2*mu_.value());

            // If the isotropic linear modulus is non-zero
            if (mag(Hp_) > SMALL)
            {
                DLambda /= 1.0 + Hp_/(3*mu_.value());

                // Update increment of yield stress
                DSigmaY = DLambda*Hp_;

                // Update yield stress
                sigmaY = sigmaYOld + DSigmaY;
            }
        }
    }
}


Foam::scalar Foam::linearElasticMisesPlasticIndep::curYieldStress
(
    const scalar curEpsilonPEq    // Current equivalent plastic strain
) const
{
    return stressPlasticStrainSeries_(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::linearElasticMisesPlasticIndep::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar            // Scaled shear modulus
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
           *curYieldStress
            (
                epsilonPEqOld + sqrtTwoOverThree_*DLambda
            );
}


void Foam::linearElasticMisesPlasticIndep::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar);

        if (i == MaxNewtonIter_)
        {
            WarningIn("linearElasticMisesPlasticIndep::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticMisesPlasticIndep::linearElasticMisesPlasticIndep
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    lambda_("zero", dimPressure, 0.0),
    stressPlasticStrainSeries_(dict),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPEq_
    (
        IOobject
        (
            "DEpsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old-time fields
    epsilon_.oldTime();
    epsilonP_.oldTime();
    epsilonPEq_.oldTime();
    plasticN_.oldTime();
    sigmaY_.oldTime();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_ = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        nu_ = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        lambda_ = E_*nu_/((1+nu_)*(1-(2*nu_)));

        // Set the bulk modulus
        K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        // Read shear modulus
        mu_ = dimensionedScalar(dict.lookup("mu"));

        // Read bulk modulus
        K_ = dimensionedScalar(dict.lookup("K"));

        // Calculate Young's modulus
        E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        // Calculate Poisson's ratio
        nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
    }
    else
    {
        FatalErrorIn
        (
            "linearElasticMisesPlasticIndep::linearElasticMisesPlasticIndep::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Check if plasticity is a nonlinear function of plastic strain
    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
    }
    else
    {
        if (stressPlasticStrainSeries_.size() == 1)
        {
            Info<< "    Perfect Plasticity" << endl;
        }
        else
        {
            Info<< "    Plasticity is linear" << endl;

            // Define linear plastic modulus
            Hp_ =
                (
                    stressPlasticStrainSeries_[1].second()
                  - stressPlasticStrainSeries_[0].second()
                )
               /(
                    stressPlasticStrainSeries_[1].first()
                  - stressPlasticStrainSeries_[0].first()
                );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticMisesPlasticIndep::~linearElasticMisesPlasticIndep()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlasticIndep::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlasticIndep::E() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            E_,
            calculatedFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlasticIndep::nu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_,
            calculatedFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlasticIndep::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            mu_,
            calculatedFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlasticIndep::lambda() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "lambda",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            lambda_,
            calculatedFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::linearElasticMisesPlasticIndep::K() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            K_,
            calculatedFvPatchScalarField::typeName
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::linearElasticMisesPlasticIndep::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField scaleFactor = 1.0 - (2.0*mu_*DLambda_/magSTrial);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::linearElasticMisesPlasticIndep::correct
(
    volSymmTensorField& sigma,
    const volTensorField& gradD,
    const volScalarField& threeKalphaDeltaT
)
{
    // Update total strain
    epsilon_ = symm(gradD);

    // Calculate deviatoric strain
    const volSymmTensorField e = dev(epsilon_);

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = 2.0*mu_*(e - dev(epsilonP_.oldTime()));

    // Calculate the yield function
    const volScalarField fTrial =
        mag(sTrial) - sqrtTwoOverThree_*sigmaY_.oldTime();

#ifdef OPENFOAMESIORFOUNDATION
    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon_.primitiveField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticN_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaY_.primitiveFieldRef();
    scalarField& DLambdaI = DLambda_.primitiveFieldRef();
    scalarField& sigmaYI = sigmaY_.primitiveFieldRef();
    const scalarField& sigmaYOldI = sigmaY_.oldTime().primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().primitiveField();
#else
    const scalar maxMagBE = max(gMax(mag(epsilon_.internalField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    scalarField& sigmaYI = sigmaY_.internalField();
    const scalarField& sigmaYOldI = sigmaY_.oldTime().internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();
#endif

    forAll(fTrialI, cellI)
    {
        // Update plasticN, DLambda, DSigmaY and sigmaY for this cell
        updatePlasticity
        (
            plasticNI[cellI],
            DLambdaI[cellI],
            DSigmaYI[cellI],
            sigmaYI[cellI],
            sigmaYOldI[cellI],
            fTrialI[cellI],
            sTrialI[cellI],
            epsilonPEqOldI[cellI],
            mu_.value(),
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
#ifdef OPENFOAMESIORFOUNDATION
        symmTensorField& plasticNP = plasticN_.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryFieldRef()[patchI];
#else
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
#endif
        const scalarField& sigmaYOldP =
            sigmaY_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            updatePlasticity
            (
                plasticNP[faceI],
                DLambdaP[faceI],
                DSigmaYP[faceI],
                sigmaYP[faceI],
                sigmaYOldP[faceI],
                fTrialP[faceI],
                sTrialP[faceI],
                epsilonPEqOldP[faceI],
                mu_.value(),
                maxMagBE
            );
        }
    }

    // Update DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;

    // Store previous iteration for residual calculation
    DEpsilonP_.storePrevIter();

    // Update DEpsilonP
    DEpsilonP_ = DLambda_*plasticN_;

    // Update total plastic strain
    epsilonP_ = epsilonP_.oldTime() + DEpsilonP_;

    // Update equivalent total plastic strain
    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq_;

    // Calculate deviatoric stress
    const volSymmTensorField s = sTrial - 2*mu_*DEpsilonP_;

    // Calculate the hydrostatic pressure
    const volScalarField p = -K_*tr(epsilon_);

    // Update the stress
    sigma = s - p*I; /*- threeKalphaDeltaT*I;*/

    // Update active yield field that is used for post-processing
    activeYield_ = pos(DEpsilonPEq_ + SMALL);
}


// ************************************************************************* //
