{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();
    
    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    
    adjustPhi(phiU, U, pd);
    
    phi = phiU +
        (
          fvc::interpolate((interface.sigmaK()+pVap)*damper)*fvc::snGrad(alphaM)
	  //+ fvc::interpolate(recoil) 
	  //+ fvc::interpolate((pVap*N)*mag(fvc::grad(alphaM))*damper)
	  - ghf*fvc::snGrad(rho)
         )*rUAf*mesh.magSf();
    
    while (pimple.correctNonOrthogonal())
    {
	pd.storePrevIter();
	
	//volScalarField rhoSource = alphaM*(fvc::ddt(rhoM)+(fvc::grad(rhoM)&U))/rhoM /*+ (1-alphaM)*fvc::ddt(rhoG)/rhoG*/;
        fvScalarMatrix pdEqn
        (
	 fvm::laplacian(rUAf, pd) == fvc::div(phi) /*+ rhoSource*/
        );
     
        pdEqn.setReference(pdRefCell, pdRefValue);
       
	pdEqn.relax();
 
        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        
        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();

	    pd.relax();

	    //U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
	    //U.correctBoundaryConditions();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
    
}
