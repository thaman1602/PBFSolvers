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

Application
    residualStressFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic thermal solid bodies. Temperature is solved and then coupled
    displacement is solved.

    Displacement field U is solved for using a total Lagrangian approach,
    also generating the strain tensor field epsilon and stress tensor
    field sigma and temperature field T.

Author
    Gowthaman Pariendhan, University College Dublin
    Philip Cardiff, University College Dublin
    Thomas Flint, The University of Manchester
   
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "thermalModel.H"
#include "FvMeshSubset.H"
#include "mapLocalToGlobal.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "linearElasticMisesPlasticIndep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "extractSolidMesh.H"
#   include "checkRegions.H"
#   include "createSubsetFields.H"
//#   include "readDivSigmaExpMethod.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
#       include "writeFields.H"

    //FatalError<< "Check Mesh implementation" << abort(FatalError);

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
#       include "readSolidMechanicsControls.H"
        int iCorr = 0;
        scalar initialResidual = 1.0;
        scalar relResD = 1.0;
        lduSolverPerformance solverPerfD;
        lduMatrix::debug = 0;

        // Solve momentum equation for displacement
        iCorr = 0;

	    //volScalarField threeKalphaDeltaT = threeKalpha*(solidResT-T0);
        threeKalphaDeltaT = threeKalpha*(solidResT-T0);
        threeKalphaDeltaT.correctBoundaryConditions();

        volVectorField gradThreeKalphaDeltaT
        (
        IOobject
        (
         "gradThreeKalphaDeltaT",
         runTime.timeName(),
         solidMesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
         ),
        solidMesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero),
        zeroGradientFvPatchScalarField::typeName
        );

        gradThreeKalphaDeltaT =
            fvc::grad(threeKalphaDeltaT, "grad(threeKalphaDeltaT)");

        surfaceVectorField threeKalphaDeltaTf =
            solidMesh.Sf()*fvc::interpolate(solidResT-T0, "deltaT");

        Info<<"impK (max/min): "<< mag(max(impK).value())<<", "<<mag(max(impK).value())<<", Dimensions: "<<impK.dimensions()<<endl;


    	/*volVectorField divThermalStresses
    	(
        IOobject
        (
         "divThermalStresses",
         runTime.timeName(),
         solidMesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
         ),
        solidMesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    	);

	divThermalStresses = fvc::div(threeKalphaDeltaT*I);*/


        Info<< "Solving for " << D.name() << nl;
        do
        {
            D.storePrevIter();
            // Linear momentum equaiton
            fvVectorMatrix DEqn
            (
                rho*fvm::d2dt2(D)
             ==
                fvm::laplacian(impKf, D, "laplacian(DD,D)")
                -fvc::laplacian(impKf, D, "laplacian(DD,D)")
                + fvc::div(sigma, "div(sigma)") - gradThreeKalphaDeltaT
                +0.1*(fvc::laplacian(impK, D, "laplacian(DD,D)") - fvc::div(impK*gradD))
            );

            /*forAll(unconstrainedRegionList, cellI)
            {
                Info<<"Reference cells: "<<unconstrainedRegionList[cellI]<<endl;
                label curCellID=unconstrainedRegionList[cellI];
            }*/

            if (unconstrainedRegionList.size())
            {
                Info<<"Fixing Reference Cells"<<endl;
                vectorField refValues (unconstrainedRegionList.size(), vector::zero);
                DEqn.setValues(unconstrainedRegionList, refValues);
            }

            solverPerfD = DEqn.solve();

            gradD = fvc::grad(D);
            mechLaw.correct(sigma, gradD,threeKalphaDeltaT);
            epsilon = symm(gradD);

            if (aitkenRelax)
            {
#               include "aitkenRelaxation.H"
            }
            else
            {
                D.relax();
            }

#           include "calculateRelResD.H"

            if (iCorr == 0)
            {
                initialResidual = solverPerfD.initialResidual();
            }

            if (iCorr % infoFrequency == 0)
            {
                Info<< "\tCorrector " << iCorr
                    << ", residual = " << solverPerfD.initialResidual()
                    << ", relative res = " << relResD;

                if (aitkenRelax)
                {
                    Info << ", aitken = " << aitkenTheta;
                }
                Info<< ", inner iters = " << solverPerfD.nIterations() << endl;
            }
        }
        while
        (
            iCorr++ == 0
         || (
                relResD > convergenceToleranceD
             && iCorr < nCorr
            )
        );

surfaceVectorField forcef = fvc::interpolate(sigma) & mesh.Sf();

	vector netForce = gSum(forcef);

//	forAll(forcef, faceI)
//	{
//	  netForce += forcef[faceI];
//	}

	// extract Values along the line x =0.05;

	const scalar currentLine = 0.05;

	DynamicList<vector> forceAlongTheLine;
	vector netForceLine = vector::zero;

	forAll(mesh.Cf(), faceI)
	{
	  if (mag(mesh.Cf()[faceI].x()-currentLine)<0.000120)
	  {
	    //Info << "coords: " << mesh.Cf()[faceI] << endl;
	    forceAlongTheLine.append(forcef[faceI]);
	    netForceLine += forcef[faceI];
	  }
	}

	Info << "Net force: " << netForce << endl;
	Info << "Net force along x=0.05: " << netForceLine << endl;

        Info<< "Solved for " << D.name()
            << " using " << solverPerfD.solverName()
            << " in " << iCorr << " iterations"
            << ", initial res = " << initialResidual
            << ", final res = " << solverPerfD.initialResidual()
            << ", final rel res = " << relResD << nl
            << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << ", ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;

//#       include "calculateEpsilonSigma.H"
#       include "writeFields.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
