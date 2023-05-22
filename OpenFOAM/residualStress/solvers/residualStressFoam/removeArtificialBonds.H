volScalarField updatedParID
(
    IOobject
    (
        "updatedParID",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    parID,
   	zeroGradientFvPatchScalarField::typeName
);

scalar nParticles = max(parID).value();

Info<<"Number of particles: "<<nParticles<<endl;

List<label> nCellsInParticles(nParticles+1, 0);

forAll(parID,cellI)
{
    if(parID[cellI]>0)
    {
        for(int i=1;i<=nParticles+1;i++)
        {
            if(parID[cellI]==i)
            {
                nCellsInParticles[i]++;
                break;
            }
        }
    }
}

//Checked upto here

forAll(updatedParID,cellI)
{
    if((solidificationTime[cellI]>0.0)||(alphaM[cellI]<0.5+SMALL))
    {
        updatedParID[cellI]=-1;
    }
}


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

List<label> unmeltedCells(nParticles+1, 0);

forAll(updatedParID,cellI)
{
    if(updatedParID[cellI]>0)
    {
        for(int i=1;i<=nParticles+1;i++)
        {
            if(updatedParID[cellI]==i)
            {
                unmeltedCells[i]++;
                unmeltedRegions[cellI]=0;
                break;
            }
        }
    }
}

updatedParID.write();
unmeltedRegions.correctBoundaryConditions();
unmeltedRegions.write();

Info<<"Total Cells in Particle 151: "<<nCellsInParticles[151]<<", unmelted Cells in particle 151: "<<unmeltedCells[151]<<endl;

//List<label> flag(nParticles+1, 1);
//volScalarField flag=1.0;

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
    scalar(1.0)
);

forAll(flag,cellI)
{
    for(int i=1;i<=nParticles+1;i++)
    {
        if(updatedParID[cellI]==i)
        {
            if(nCellsInParticles[i]==unmeltedCells[i])
            {
                flag[cellI]=0;
            }
        }
    }
}

flag.correctBoundaryConditions();
flag.write();

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

forAll(region,cellI)
{
    if(flag[cellI]==1 && alphaM[cellI]>0.5)
    {
        label neiCount = 0;
        label plateCheck = 0;
        const unallocLabelList& nei = mesh.cellCells()[cellI];
        forAll(nei,neiI)
        {
            if((alphaI[nei[neiI]]>0.5+SMALL) && (flag[nei[neiI]] == 1))
            {
                neiCount++;
            }
            if (updatedParID[nei[neiI]]==0 && updatedParID[cellI]>0)
            {
                plateCheck=1;
            }
        }
        if(neiCount>0 && (plateCheck==0 || updatedParID[cellI]==0))
        {
            region[cellI]=1; //Add solid cells that have atleast one neighbour. Floating cells cause floating point exception if all BC are traction free
            extractedRegion[cellI]=1;

        }
    }
}

extractedRegion.correctBoundaryConditions();
extractedRegion.write();


/*forAll (alphaI, cellI)
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
            }
        }
}*/



