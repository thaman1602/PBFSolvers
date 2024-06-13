volScalarField normalFlag
(
    IOobject
    (
        "normalFlag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

label counter = 1;
DynamicList<label> interface;
forAll (solidificationTime, cellI)
{
    if(alphaM[cellI] > 0.5 && solidificationTime[cellI] > SMALL)
    {
        label neiCount = 0;
        const unallocLabelList& nei = mesh.cellCells()[cellI];
        forAll(nei,neiI)
        {
            if(alphaM[nei[neiI]]>0.5 && solidificationTime[nei[neiI]]==-1)
            {
                neiCount++;
            }
        }
        if(neiCount>0)
        {
            normalFlag[cellI]=counter;
            counter++;
            interface.append(cellI);
        }
    }
    if(alphaM[cellI]>0.5 && solidificationTime[cellI]<SMALL)
    {
        normalFlag[cellI]=-1;
    }
}

counter = interface.size();
while(!interface.empty())
{
    label curCellId = interface.remove(); //front.removeHead();
    const labelList& faces = cells[curCellId];
        
    forAll(faces, faceI)
    {
        const label& face = faces[faceI];
        label neighbourID = owner[face] == curCellId ? neighbour[face] : owner[face];

        if(solidificationTime[neighbourID]<SMALL && normalFlag[neighbourID]<SMALL)
        {  
            scalar deltaT = mag((C[curCellId] - C[neighbourID]) & gradTSolf[face]);
              
            resT[neighbourID] = resT[curCellId] - deltaT;
            normalFlag[neighbourID]=normalFlag[curCellId]+1000;
        }
    }
    counter--;
}

normalFlag.write();