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
    porosityCalculator

Description
    Postprocessing utility to compute the porosity in the substrate.

Author
    Gowthaman Parivendhan, University College Dublin

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    
    instantList Times = runTime.times();
    
    if (Times.size())
    {
        label startTime = Times.size()-1;
        runTime.setTime(Times[startTime], startTime);
    }
    else
    {
        runTime.setTime(instant(0, runTime.constant()), 0);
    }

    Info<< "Time = " << runTime.timeName() << endl;

    IOdictionary laserProperties
    (
        IOobject
        (
            "laserProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word laserPatchName = word(laserProperties.lookup("laserPatchName"));

    volScalarField flag
        (
        IOobject
        (
            "flag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0",dimless, 1.0),
        zeroGradientFvPatchScalarField::typeName
        );


    volScalarField alphaM
    (
        IOobject
        (
            "alpha.material",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField solidificationTime
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

    Info<<"Creating cell columns."<<endl;

#   include "createCellColumns.H"

    Info<<"Tracking Top surface of the powder bed."<<endl;
    forAll(cellColumns, colI)
    {
        DynamicList<label>& curCol = cellColumns[colI];
        scalar rayFraction = 1.0;
    
        forAll(curCol,colJ)
        {
            label cellID = curCol[colJ];
            if(pos(rayFraction))
            {
                if(alphaM[cellID]>(0.95+SMALL))
                {
                    rayFraction = -1.0;
                }
                flag[cellID]=0.0;
            }
        }
    }

    scalar totalGasVolume = 0.0;
    scalar totalMeltTrackVolume = 0.0;

    forAll(alphaM,cellI)
    {
        scalar volCell = mesh.V()[cellI];
        scalar gasVolume = 0.0;
        scalar cellVolume = 0.0;
        if(solidificationTime[cellI]>0.0)    
        {
            gasVolume = flag[cellI] * (1-alphaM[cellI]) * volCell;
            cellVolume = flag[cellI] * volCell;
        }

        totalGasVolume += gasVolume;
        totalMeltTrackVolume += cellVolume;
    }

    scalar porosity=(totalGasVolume/totalMeltTrackVolume)*100.0;

    Info<<"Porosity of the substrate = "<<porosity<<"%"<<endl;

    Info<<"End."<<endl;
    return 0;
}


// ************************************************************************* //
