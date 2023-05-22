    Info<< "Reading field D\n" << endl;
    volVectorField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    volTensorField gradD
    (
        IOobject
        (
            "grad(D)",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    surfaceVectorField snGradD
    (
        IOobject
        (
            "snGrad(D)",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            solidMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        solidMesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    Info<< "Reading field resT\n" << endl;
    volScalarField solidResT
    (
        IOobject
        (
            "solidResT",
            runTime.timeName(),
            solidMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        solidMesh,
        dimensionedScalar("zero",  dimensionSet(0,0,0,1,0,0,0), 0.0),
    	zeroGradientFvPatchScalarField::typeName
    );

    solidResT = subsetMesh.interpolate(resT);
    solidResT.correctBoundaryConditions();

    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            solidMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        solidMesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    Info<< "Reading thermal model\n" << endl;
    thermalModel thermal(solidResT);

    surfaceVectorField n = solidMesh.Sf()/solidMesh.magSf();

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

    linearElasticMisesPlasticIndep mechLaw(solidMesh, rheologyProperties);

    volScalarField rho = mechLaw.rho();

    volScalarField impK
    (
        IOobject
        (
            "impK",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidMesh,
        dimensionedScalar("zero",  dimensionSet(1,-1,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    impK=mechLaw.impK();

    surfaceScalarField impKf = fvc::interpolate(impK, "impK");

    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidMesh,
        dimensionedScalar("zero",  dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    mu=mechLaw.mu();

    volScalarField lambda
    (
        IOobject
        (
            "lambda",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidMesh,
        dimensionedScalar("zero",  dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    lambda=mechLaw.lambda();

    volScalarField threeKalpha
    (
        IOobject
        (
            "threeKalpha",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidMesh,
        dimensionedScalar("zero",  dimensionSet(1,-1,-2,-1,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    threeKalpha=scalar(3.0)*mechLaw.K()*thermal.alpha();
    threeKalpha.correctBoundaryConditions();

    surfaceScalarField threeKalphaf =
        fvc::interpolate(threeKalpha, "threeKalpha");

    volScalarField T0 = thermal.T0();

    volScalarField threeKalphaDeltaT
    (
        IOobject
        (
            "threeKalphaDeltaT",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidMesh,
        dimensionedScalar("zero",  dimensionSet(1,-1,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField k
    (
        "DT",
        thermal.k()
     ); //Dummy variable for fvPatchFields - set to 0.0 

    // for aitken relaxation
    volVectorField aitkenDelta
    (
        IOobject
        (
         "aitkenDelta",
         runTime.timeName(),
         solidMesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
         ),
        solidMesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // aitken relaxation factor
    scalar aitkenInitialRes = 1.0;
    scalar aitkenTheta = 0.01;
