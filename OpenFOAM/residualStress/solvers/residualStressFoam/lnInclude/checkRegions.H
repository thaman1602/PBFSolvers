regionSplit rs(solidMesh);
DynamicList<label> unconstrainedRegionList(solidMesh.nCells()); //List to store the free regions - eg: group of cells separated from the base region 
//DynamicList<label> floatingCellsList(rs.nRegions()); //Single floating cells that are to be removed from the mesh 

if (rs.nRegions() == 1)
{
    Info<< "    Number of regions: " << rs.nRegions() << " (OK)."
        << endl;

}
else
{
    Info<< "   *Number of regions: " << rs.nRegions() << endl;

    Info<< "    The mesh has multiple regions which are not connected "
    "by any face." << endl;

    // Count number of cells in all regions
    labelList nCellsInRegions(rs.nRegions(), 0);
    PtrList<DynamicList<label>> cellsInRegions(rs.nRegions());

    forAll (rs, rsI)
    {
        nCellsInRegions[rs[rsI]]++;
    }

    forAll (cellsInRegions, regI)
    {
        cellsInRegions.set(regI, new DynamicList<label>(nCellsInRegions[regI]));
    }

    forAll (rs, rsI)
    {
        DynamicList<label>& curReg = cellsInRegions[rs[rsI]];
        curReg.append(rsI);
    }

    forAll (cellsInRegions, regI)
    {
        DynamicList<label>& curReg = cellsInRegions[regI];
        Info<<"Region: "<<regI<<", number of cells: "<<curReg.size()<<endl;
        forAll(curReg,cellI)
        {
            /*if(curReg.size()==1)
            {
                floatingCellsList.append(curReg[cellI]);
            }
            else
            {*/
                unconstrainedRegionList.append(curReg[cellI]); //Take reference cell from each region to fix geometry
                break;
            //}
        }
    }        

    Info<< "Number of cells per region: " << nl;

    forAll (nCellsInRegions, regionI)
    {
        Info<< tab << regionI << tab << nCellsInRegions[regionI] << nl;
    }

    Info<< endl;
}

unconstrainedRegionList.shrink();
//floatingCellsList.shrink();

Info<<"unconstrainedRegionList size: "<<unconstrainedRegionList.size()<<endl;
//Info<<"floatingCellsList size: "<<floatingCellsList.size()<<endl;

/*removeCells cellRemover(solidMesh);

Info<< "Mesh has:" << solidMesh.nCells() << " cells."
    << " Removing:" << floatingCellsList.size() << " cells" << endl;

labelList exposedFaces(cellRemover.getExposedFaces(floatingCellsList));

directTopoChange meshMod(solidMesh);
cellRemover.setRefinement
(
        floatingCellsList,
        exposedFaces,
        labelList(exposedFaces.size(), patchI),
        meshMod
);

Info<< "Morphing ..." << endl;

const Time& runTime = solidMesh.time();

autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(solidMesh, false);
cellRemover.updateMesh(morphMap());*/

forAll(unconstrainedRegionList,cellI)
{
    Info<<"Cell ID: "<<unconstrainedRegionList[cellI]<<endl;

}



