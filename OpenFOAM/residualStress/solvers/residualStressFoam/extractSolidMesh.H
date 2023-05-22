fvMeshSubset subsetMesh
(
    IOobject
    (
        "solid",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

labelList region(mesh.nCells(), 0);

volScalarField alphaM
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

volScalarField resT
(
    IOobject
    (
        "resT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
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

volScalarField unmeltedRegions
(
    IOobject
    (
        "unmeltedRegions",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField extractedRegion
(
    IOobject
    (
        "extractedRegion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(0.0),
    zeroGradientFvPatchScalarField::typeName
);


dictionary simType = mesh.solutionDict().subDict("simType");

bool PBF = simType.lookupOrDefault<bool>("PBF", false);

scalarField& alphaI = alphaM.internalField();

//checkInterfacePatchName(mesh);
label solidInterfaceID = mesh.boundaryMesh().findPatchID("solidInterface");

Info<<"solidInterfaceID: "<<solidInterfaceID<<endl;


if(PBF)
{
    volScalarField parID
    (
        IOobject
        (
            "parID",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    #include "removeArtificialBonds.H"
}
else
{
    forAll (alphaI, cellI)
    {
        if(alphaI[cellI] > 0.5)
        {
	       label neiCount = 0;
	       const unallocLabelList& nei = mesh.cellCells()[cellI];
	       forAll(nei,neiI)
	       {
                if(alphaI[nei[neiI]]>0.5)
                {
                    neiCount++;
                }
	       }
            if(neiCount>0)
            {
                region[cellI]=1; //Add solid cells that have atleast one neighbour. Floating cells cause floating point exception
                extractedRegion[cellI]=1; 
            }
           /*if(solidificationTime[cellI]>0+SMALL)
            {
                region[cellI]=1; //Add solid cells that have atleast one neighbour. Floating cells cause floating point exception
                extractedRegion[cellI]=1; 
            }*/
        }
    }  
}
subsetMesh.setLargeCellSubset(region, 1, solidInterfaceID);

fvMesh& solidMesh = subsetMesh.subMesh();

