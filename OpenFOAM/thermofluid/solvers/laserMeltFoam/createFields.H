    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha\n" << endl;
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

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
	IOobject
	(
	    "T",
	     runTime.timeName(),
	     mesh,
	     IOobject::MUST_READ,
	     IOobject::AUTO_WRITE
	),
	mesh
    );

#   include "readThermalProperties.H"


    volScalarField cpS
    (
    IOobject
    (
        "cpS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(0,2,-2,-1,0,0,0), 0.0),
    zeroGradientFvPatchScalarField::typeName
    );

    volScalarField gT
    (
        IOobject
        (
            "gT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        liquidFraction(T, Tl, Ts),
        zeroGradientFvPatchScalarField::typeName
     );

    Info << "After gT update:\t" << "  Min(gT) = " << min(gT).value()
    << "  Max(gT) = " << max(gT).value() << endl;

    Info<< "Reading field alpha.liquid\n" << endl;
    volScalarField alphaL
    (
        IOobject
        (
            "alpha.liquid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        alphaM * gT,
        zeroGradientFvPatchScalarField::typeName
     );

    Info << "After alphaL update:\t" << "  Min(alphaL) = " << min(alphaL).value()
    << "  Max(alphaL) = " << max(alphaL).value() << endl;

#   include "createPhi.H"
	
    Info<< "Reading transportProperties\n" << endl;
    threePhaseMixture threePhaseProperties(U, phi, "alpha.material", "gT","T");

    volScalarField rhoM
    (
     IOobject
     (
      "rhoM",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar("0", dimensionSet(1,-3,0,0,0,0,0), 7950.0),
     zeroGradientFvPatchScalarField::typeName
     );

    const dimensionedScalar& rhoL = threePhaseProperties.rhoL();
    const dimensionedScalar& rhoS = threePhaseProperties.rhoS();
    const dimensionedScalar& rhoG = threePhaseProperties.rhoG();

    volScalarField cM
    (
        IOobject
        (
        "cM",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT
         ),
     mesh,
     dimensionedScalar("0", dimensionSet(0,2,-2,-1,0,0,0), 0.0),
     zeroGradientFvPatchScalarField::typeName
     );

    volScalarField kM
    (
        IOobject
        (
         "kM",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1,1,-3,-1,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField kEff
    (   
        IOobject
        (
        "kEff",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT
         ), 
     alphaM*(kM)+(scalar(1)-alphaM)*kG,
     zeroGradientFvPatchScalarField::typeName
     );

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
	    IOobject::NO_WRITE
         ),
        alphaM*((gT*rhoL)+((scalar(1)-gT)*rhoS)) + (scalar(1) - alphaM)*rhoG,
        zeroGradientFvPatchScalarField::typeName
     );
    rho.oldTime();

    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
         ),
        rho*(alphaM*cM + (scalar(1)-alphaM)*cpG),
        zeroGradientFvPatchScalarField::typeName
     );

    volScalarField A
    (
     IOobject
     (
      "A",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT
      ),
     -1.0 * C * alphaM * (Foam::pow(scalar(1)-gT,2))/(Foam::pow(gT,3)+SMALL),
     zeroGradientFvPatchScalarField::typeName
     );

    volScalarField rhoT
    (
        IOobject
        (
            "rhoT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
        mesh,
        dimensionedScalar("0", dimensionSet(1,-3,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
     );

    #include "updateProperties.H"

    Info<<"max/min(rho): "<<max(rho).value()<<", "<<min(rho).value()<<endl;

    Info<<"max/min(cpS): "<<max(cpS).value()<<", "<<min(cpS).value()<<endl;
    Info<<"max/min(cM): "<<max(cM).value()<<", "<<min(cM).value()<<endl;

    Info<<"max/min(kM): "<<max(kM).value()<<", "<<min(kM).value()<<endl;
    Info<<"max/min(kEff): "<<max(kEff).value()<<", "<<min(kEff).value()<<endl;

    volScalarField deltaH
    (
    IOobject
    (
    "deltaH",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ),
    rhoM*(cpL - cpS)*T + rhoM*Hf
    );

    Info<< "Creating field implicit source Sp\n" << endl;
    volScalarField SI
    (
     IOobject
     (
      "SI",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar("0", dimensionSet(1,-1,-3,-1,0,0,0), 0.0),
     zeroGradientFvPatchScalarField::typeName
    );

    Info<< "Creating field explicit source Sc\n" << endl;
    volScalarField SE
    (
     IOobject
     (
      "SE",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar("0", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
     zeroGradientFvPatchScalarField::typeName
    );

    dictionary controls = mesh.solutionDict().subDict("Controls");

    double relax =
    controls.lookupOrDefault<double>("relaxCorr", 1.0);

    label nTCorrectors(readLabel(pimple.dict().lookup("nTCorrectors")));
    scalar temperatureTol
    (
     readScalar(pimple.dict().subDict("residualControl").
            subDict("temperatureTol").lookup("tolerance"))
    );

    scalar alternativeTol
    (
     readScalar(pimple.dict().subDict("residualControl").
            subDict("alternativeTol").lookup("tolerance"))
    );

    scalar absTol
    (
     readScalar(pimple.dict().subDict("residualControl").
            subDict("absTol").lookup("tolerance"))
     );

    scalar gTTol
    (
     readScalar(pimple.dict().subDict("residualControl").
            subDict("gTTol").lookup("tolerance"))
     );

    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alphaM solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
     fvc::interpolate(rho)*phi
    );

    surfaceScalarField rhoCpPhi
    (
        IOobject
        (
            "rhoCpPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
        fvc::interpolate(rhoCp)*phi
    );

    surfaceScalarField rhoMPhi
    (
        IOobject
        (
            "rhoMPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
     fvc::interpolate(alphaM*((gT*rhoL)+((scalar(1)-gT)*rhoS)))*phi*Hf
     );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + rho*gh
    );

    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, pimple.dict(), pdRefCell, pdRefValue);
    mesh.schemesDict().setFluxRequired(pd.name());

    scalar pRefValue = 0.0;

    if (pd.needReference())
    {
        pRefValue = readScalar(pimple.dict().lookup("pRefValue"));

        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pdRefCell)
        );
    }

    // Construct interface from alphaM distribution
    interfaceProperties interface(alphaM, U, threePhaseProperties);

    const dimensionedScalar DSigbyDT (threePhaseProperties.lookup("STgrad")); // Surface Tension gradient

     // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, threePhaseProperties)
    );


    volScalarField pVap
    (
        IOobject
        (
            "pVap",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
        mesh,
        dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
     ); // Recoil Pressure

    volScalarField Qv
    (
        IOobject
        (
            "Qv",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
        mesh,
        dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
     );
     // Evaporation Cooling

    //Residual Stress calculation

    volVectorField gradTSol
    (
        IOobject
        (
	 "gradTSol",
	 runTime.timeName(),
	 mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
	 ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
     );

    volScalarField solidificationTime
    (
        IOobject
        (
	 "solidificationTime",
	 runTime.timeName(),
	 mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
	 ),
        mesh,
        dimensionedScalar("minusOne", dimless, -1),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField ddtGT = fvc::ddt(gT);
