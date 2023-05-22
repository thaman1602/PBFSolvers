if (runTime.outputTime())
{
    volScalarField epsilonEq
    (
        IOobject
        (
            "epsilonEq",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    );
    
    volScalarField globalEpsilonEq
    (
        IOobject
        (
            "globalEpsilonEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",  dimless, 0.0)
     );
    
    mapLocalToGlobal(globalEpsilonEq,epsilonEq,subsetMesh);
    
    Info<< "Max epsilonEq = " << max(epsilonEq).value() << endl;

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    volScalarField globalSigmaEq
    (
     IOobject
     (
      "globalSigmaEq",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
      mesh,
      dimensionedScalar("zero", dimForce/dimArea, 0.0)
     );

    mapLocalToGlobal(globalSigmaEq,sigmaEq,subsetMesh);
  
    Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;

    volVectorField globalD
    (
     IOobject
     (
      "globalD",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
      mesh,
     dimensionedVector("zero", dimless, vector::zero)
    );
    
    volSymmTensorField globalEpsilon
    (
     IOobject
     (
      "globalEpsilon",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedSymmTensor("zero", dimless, symmTensor::zero)
     );
    
    volSymmTensorField globalSigma
    (
     IOobject
     (
      "globalSigma",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     mesh,
     dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
     );
    
    mapLocalToGlobal(globalD,D,subsetMesh);
    mapLocalToGlobal(globalEpsilon,epsilon,subsetMesh);
    mapLocalToGlobal(globalSigma,sigma,subsetMesh);
    
    runTime.write();
    extractedRegion.write();
    unmeltedRegions.write();
}
