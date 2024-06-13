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
        dimensionedScalar("zero", dimless, -1),
        zeroGradientFvPatchScalarField::typeName
    );

    vector solNormal(0,0,0);

    label counter = 1;
    DynamicList<label> interface;
    DynamicList<label> solidList;

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
                if(solidificationTime[cellI]>SMALL)
                {
                    interface.append(cellI);
                }
           }
        }
        if(alphaM[cellI]>0.5 && solidificationTime[cellI]<SMALL)
        {
            solidList.append(cellI);
            normalFlag[cellI]=-1;
        }
    } 

    scalar minSolidificationTime = GREAT;
    /*scalar coldCell = -1;
    forAll(solidificationTime, cellI)
    {
        if(solidificationTime[cellI]>0)
        {
            if(solidificationTime[cellI]<minSolidificationTime)
            {
                minSolidificationTime=solidificationTime[cellI];
                coldCell = cellI;
            }
        }
    }
    Info<<"Cold cell solidification time: "<<minSolidificationTime<<endl;
    interface.append(coldCell);*/

    List<vector> firstCells(interface.size(),vector::zero);
    List<vector> lastCells(interface.size(),vector::zero);
    List<scalar> interfaceResT(interface.size(),0.0);
    counter = interface.size();
    //normalFlag[coldCell]=counter;
    //vector firstCell = vector::zero;
    //vector lastCell = vector::zero;


    while(!interface.empty())
    {
        label curCellID = interface.remove(); 
        vector solNormal(0,0,0);

        solNormal=gradTSol[curCellID]/mag(gradTSol[curCellID]);
        vector curCellLoc = mesh.C()[curCellID];
        //firstCell = curCellLoc;
        firstCells[counter-1]=curCellLoc;
        interfaceResT[counter-1]=resT[curCellID];
        vector forwardNewPoint = vector::zero;
        do
        {
            //solNormal=gradTSol[curCellID]/mag(gradTSol[curCellID]);
            forwardNewPoint = curCellLoc - (1e-5+SMALL)*solNormal;
            curCellID = mesh.findCell(forwardNewPoint);
            //normalFlag[curCellID] = 1.0;
            normalFlag[curCellID] = counter;
            resT[curCellID] = interfaceResT[counter-1];
            curCellLoc=forwardNewPoint;
            lastCells[counter-1]=curCellLoc;
        }
        while(curCellLoc.component(0)<2e-3 && curCellLoc.component(2)<2e-3 && solidificationTime[curCellID]<0);
        counter--;
    }

    //normalFlag.write();
    label oldFlag = 0;
    scalar oldResT = -467.588;
    for (int i = 0; i < firstCells.size(); ++i)
    {
        vector p1 = firstCells[i];
        vector p2 = lastCells[i];
        //vector p1 = firstCell;
        //vector p2 = lastCell;
        //Info<<"Interface ResT List: "<<interfaceResT[i]<<endl;
        scalar slope = (p1.component(2)-p2.component(2))/(p1.component(0)-p2.component(0));

        //while(!solidList.empty())
        //{
        forAll(solidList,listI)
        {
            scalar trialCoord = 0.0;
            label curCellID = solidList[listI];
            vector curCellLoc = mesh.C()[curCellID];
            trialCoord=slope*(curCellLoc.component(0) - p1.component(0))+p1.component(2);

            if(trialCoord<curCellLoc.component(2) && normalFlag[curCellID]<0)
            {
                //if(resT[curCellID]<-GREAT)
                //{
                //normalFlag[curCellID]=(oldFlag+flags[i])/2;
                resT[curCellID]=(oldResT+interfaceResT[i])/2;
                //solidList.remove();
            }
            //}
        }
        //oldFlag=flags[i];
        oldResT=interfaceResT[i];
    }
        
    forAll(solidList,listI)
    {
        label curCellID = solidList[listI];
        if(resT[curCellID]<-GREAT /*normalFlag[curCellID]<0*/)
        {
            //normalFlag[curCellID]=oldFlag;
            resT[curCellID]=oldResT;    
        }   
    }

normalFlag.write();