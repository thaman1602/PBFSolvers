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
    TResIntegrator

Description
    Utility to compute the TResIntegrator.

    Procedure based on Schilling et al, Approach on simulation of
    solidification and shrinkage of gravity cast salt cores, Simulation
    Modelling Practice and Theory 107 (2021) 102231.

Author
    Gowthaman Parivendhan, University College Dublin
    Philip Cardiff, University College Dublin
    Tom Flint, The University of Manchester

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sortDynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList times = runTime.times();

    #include "createMesh.H"

    if (times.size())
    {
        const label startTime = times.size() - 1;
        runTime.setTime(times[startTime], startTime);
    }
    else
    {
        runTime.setTime(instant(0, runTime.constant()), 0);
    }

    volScalarField resT
    (
        IOobject
        (
            "resT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("VGREAT", dimTemperature, -VGREAT),
        zeroGradientFvPatchScalarField::typeName
    );


    Info<< "Time = " << runTime.timeName() << endl;

    const volVectorField gradTSol
    (
        IOobject
        (
            "gradTSol",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const surfaceVectorField gradTSolf(fvc::interpolate(gradTSol));

    const volScalarField solidificationTime
    (
        IOobject
        (
            "solidificationTime",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField meltRegions
    (
        IOobject
        (
            "meltRegions",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );

    const volScalarField alphaM
    (
        IOobject
        (
            "alpha.material",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const scalarField& solidificationTimeI = solidificationTime.internalField();

    Info<<"solidificationTime (min/max)"<< max(solidificationTime).value()
        <<", "<<min(solidificationTime).value()<<endl;

    const volVectorField& C = mesh.C();
    const cellList& cells = mesh.cells();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    label cellCounter = 1;
    forAll(solidificationTime, cellI)
    {
        if (solidificationTime[cellI] < SMALL)
        {
            meltRegions[cellI] = -1;
            resT[cellI] = 0.0;
            cellCounter++;
        }
    }

    Info<< "Cell counter: " << cellCounter << endl;

    DynamicList<label> front;

    while (gMax(meltRegions) > -1)
    {
        dimensionedScalar maxSolidificationTime = -GREAT;

        label hotCellID = -1;
        forAll(cells,cellI)
        {
            if (meltRegions[cellI] > -1)
            {
                if (solidificationTime[cellI] > maxSolidificationTime.value())
                {
                    maxSolidificationTime.value() = solidificationTime[cellI];
                    hotCellID = cellI;
                }
            }
        }

        Info<< "maxSolidificationTime: " << maxSolidificationTime << endl;

        front.append(hotCellID);

        Info<< "Hottest Cell: " << hotCellID << endl;

        resT[hotCellID] = 0;

        Info<<"SolTime hot cell: " << solidificationTimeI[hotCellID] << endl;

        while (!front.empty())
        {
            label cellIMax = front.remove();
            meltRegions[cellIMax] = -1;
            const labelList& faces = cells[cellIMax];

            forAll(faces, faceI)
            {
                const label& face = faces[faceI];
                if (face < mesh.nInternalFaces())
                {
                    const label cellI =
                        owner[face] == cellIMax ? neighbour[face] : owner[face];

                    if
                    (
                        resT[cellI] < -GREAT
                     && solidificationTime[cellI] > SMALL
                    )
                    {
                        const scalar deltaT =
                            mag
                            (
                                (C[cellIMax] - C[cellI]) & gradTSolf[face]
                            );

                        if
                        (
                            solidificationTime[cellIMax]
                         >= solidificationTime[cellI]
                        )
                        {
                            resT[cellI] = resT[cellIMax] - deltaT;
                        }
                        else
                        {
                            resT[cellI] = resT[cellIMax] + deltaT;
                        }
                        front.append(cellI);

                        if (front.size() > 1)
                        {
                            sortDynamicList(front, solidificationTime);
                        }
                    }
                }
            }
        }
    }

    Info<< "resT (min/max)" << min(resT).value() << ", " << max(resT).value() <<endl;

    forAll(mesh.C(), cellI)
    {
        if
        (
            (solidificationTime[cellI] < SMALL)
         && (resT[cellI] < -GREAT)
        )
        {
            resT[cellI] = 0.0;
        }
    }

    resT.correctBoundaryConditions();

    const scalar minResT = mag(min(resT).value());

    forAll(mesh.C(), cellI)
    {
        if (solidificationTime[cellI] > SMALL || resT[cellI] < 0)
        {
            resT[cellI] = -1.0*(resT[cellI] + minResT);
        }
    }

    resT.correctBoundaryConditions();

    Info<< "resT (min/max)" << min(resT).value() << ", " << max(resT).value()
        <<endl;

    resT.write();

    Info<< "End" << endl << endl;

    return 0;
}


// ************************************************************************* //
