const dictionary& stressControl =
  mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));

scalar convergenceToleranceD(readScalar(stressControl.lookup("D")));

int infoFrequency(readInt(stressControl.lookup("infoFrequency")));

Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
