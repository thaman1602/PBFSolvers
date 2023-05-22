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

Class
    linearElasticMisesPlasticIndep

Description
    Elasto-plastic mechanical law as described by Simo & Hughes (1998)
    in Box 3.2.

    Linear (Hookean) elasticity with Mises/J2 plasticity.

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    More details found in:

    Simo & Hughes, Computational Inelasticity, 1998, Springer.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

    This law has been extracted from solids4foam so that it can be used
    independently.

SourceFiles
    linearElasticMisesPlasticIndep.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMisesPlasticIndep_H
#define linearElasticMisesPlasticIndep_H

#include "volFields.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMisesPlasticIndep
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Initial density
        const dimensionedScalar rho_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        dimensionedScalar lambda_;

        //- Table of post-yield stress versus plastic strain
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Yield stress
        volScalarField sigmaY_;

        //- Incremental change in sigmaY
        volScalarField DSigmaY_;

        //- Total strain
        volSymmTensorField epsilon_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;

        //- Equivalent plastic strain
        volScalarField epsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        //- An iterative procedure is used when the yield stress is a nonlinear
        //  function of plastic strain
        const bool nonLinearPlasticity_;

        //- Linear plastic modulus. It is only used if plasticity is linear,
        //  defined by two points on the stress-plastic strain curve
        scalar Hp_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticMisesPlasticIndep(const linearElasticMisesPlasticIndep&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticMisesPlasticIndep&);

        //- Update plasticN, DLambda, DSigmaY and sigmaY
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Update plasticN, DLambda, DSigmaY and sigmaY
        void updatePlasticity
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
            const scalar maxMagDEpsilon    // Max strain increment magnitude
        ) const;

        //- Return the current yield stress
        scalar curYieldStress
        (
            const scalar curEpsilonPEq    // Current equivalent plastic strain
        ) const;

        //- Evaulate current value of the yield function
        scalar yieldFunction
        (
            const scalar epsilonPEqOld,    // Old equivalent plastic strain
            const scalar magSTrial,        // Deviatoric trial stress magnitude
            const scalar DLambda,          // Plastic multiplier
            const scalar muBar            // Scaled shear modulus
        ) const;

        //- Iteratively calculate plastic multiplier increment (DLambda)
        //  and current yield stress using Newton's method
        void newtonLoop
        (
            scalar& DLambda,               // Plastic multiplier
            scalar& curSigmaY,             // Current yield stress
            const scalar epsilonPEqOld,    // Old equivalent plastic strain
            const scalar magSTrial,        // Deviatoric trial stress magnitude
            const scalar muBar,            // Scaled shear modulus
            const scalar maxMagDEpsilon    // Max strain increment magnitude
        ) const;


public:

    //- Runtime type information
    TypeName("linearElasticMisesPlasticIndep");

    // Constructors

        //- Construct from dictionary
        linearElasticMisesPlasticIndep
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~linearElasticMisesPlasticIndep();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        virtual tmp<volScalarField> E() const;

        virtual tmp<volScalarField> nu() const;

        virtual tmp<volScalarField> mu() const;

        virtual tmp<volScalarField> lambda() const;

        //-Return Bulk modulus
        virtual tmp<volScalarField> K() const;   

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Update the stress
        virtual void correct
        (
            volSymmTensorField& sigma,
            const volTensorField& gradD,
            const volScalarField& gradThreeKalphaDeltaT
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
