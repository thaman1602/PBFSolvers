IOdictionary rheologyProperties
(
 IOobject
 (
  "rheologyProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );

const scalar ht(readScalar(rheologyProperties.lookup("ht")));

volScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
); //Heaviside function



label counter = 1;
DynamicList<vector> interfaceList;
forAll (solidificationTime, cellI)
{
    if(alphaM[cellI] > 0.5 && solidificationTime[cellI] > SMALL)
    {
        const labelList& faces = cells[cellI];
        
        forAll(faces, faceI)
        {
            const label& face = faces[faceI];
            label faceOwner = owner[face];
            label faceNeighbour = neighbour[face];
            if(face < mesh.nInternalFaces())
            {
            if(alphaM[faceOwner]>0.5 && alphaM[faceNeighbour]>0.5)
            {
                if(solidificationTime[faceOwner] == -1 || solidificationTime[faceNeighbour]==-1)
                {
                    Info<< counter<<": "<<mesh.Cf()[face] << endl;
                    interfaceList.append(mesh.Cf()[face]);
                    counter++;
                }
            }
            }
        }
        /*label neiCount = 0;
        const unallocLabelList& nei = mesh.cellCells()[cellI];
        forAll(nei,neiI)
        {
            if(alphaM[nei[neiI]]>0.5 && solidificationTime[nei[neiI]]==-1)
            {
                neiCount++;
                resT[nei[neiI]]=resT[cellI]/1.5;
            }
        }*/
        /*if(neiCount>0)
        {
            normalFlag[cellI]=counter;
            counter++;
            interface.append(cellI);
        }*/
    }
    /*if(alphaM[cellI]>0.5 && solidificationTime[cellI]<SMALL)
    {
        normalFlag[cellI]=-1;
    }*/
}

scalar distance = GREAT;
forAll(phi,cellI)
{
    forAll(interfaceList, listI)
    {
        phi[cellI]=mag(interfaceList[listI]-C[cellI]);
        if(phi[cellI]<distance)
        {
            distance=phi[cellI];
        }
    }
    phi[cellI]=distance;
    distance=GREAT;
}

forAll(H,cellI)
{
    if(phi[cellI]<=ht && phi[cellI]>=(-1*ht))
    {
        H[cellI]=(0.5*phi[cellI]/ht)+(0.5*Foam::sin(M_PI*phi[cellI]/ht)/M_PI);
    }
    if(phi[cellI]>ht)
    {
        H[cellI]=1;
    }
    if(phi[cellI]<(-1*ht))
    {
        H[cellI]=0;
    }
}

H.write();
phi.write();
