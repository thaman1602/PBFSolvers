/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    laserMeltFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids and using a VOF
    (volume of fluid) phase-fraction based interface capturing approach and solidification of liquid based on Voller method.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e.  laminar, RAS or LES may be selected.
 
Author
    Gowthaman Parivendhan, UCD
 
        "When I wrote this, only God and I understood what I was doing
        Now, only God knows"
                                                    -Karl Weierstrass
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "isoAdvection.H"
#include "laserThreePhaseMixture.H"
#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "interpolationTable.H"
#include "liquidFractionFunctions.H"
#include "bound.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "readLaserProperties.H"
#   include "createCellColumns.H"
#   include "createTimeControls.H"
#   include "correctPhi.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
    	#include "readTimeControls.H"
    	#include "CourantNo.H"
        #include "alphaCourantNo.H"
        
        // Set the timestep based on the maximum alpha courant / courant number
        #include "setAlphaDeltaT.H" // From atomizationDNS implementation - Vuko Vukcevic

        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

	//#include "updateLaser.H"
 
        // Pressure-velocity corrector
        while (pimple.loop())
        {
            threePhaseProperties.correct(); // User library to read and calculate three phase properties, values provided in transportProperties Directory
	    
	    ddtGT.oldTime();

            // Advection equation and update the coefficients
            #include "alphaEqn.H"
	 
	    #include "updateLaser.H"
	    //#include "TEqn.H"

	    //#include "updateProperties.H"	

            //Set UEqn matrix, momentum predictor turned off in all solidification cases
            #include "UEqn.H"

            
	    //#include "TEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

	    #include "TEqn.H"

            #include "continuityErrs.H"
	    
            p = pd + rho*gh;

            if (pd.needReference())
            {
                p += dimensionedScalar
                (
                    "p",
                    p.dimensions(),
                    pRefValue - getRefCellValue(p, pdRefCell)
                );
                pd = p - rho*gh;
            }
            
	    #include "updateProperties.H"
	    
            turbulence->correct();
            Info << "\n" << endl;
        }

	volVectorField gradT = fvc::grad(T);

	forAll(gT, cellI)
	{
	    if(T[cellI]<Ts.value() && ddtGT.oldTime()[cellI]<0.0 && gT.oldTime()[cellI]>(0.0+SMALL))
	      {
	          scalar t_sol = GREAT;
		  if (mag(ddtGT.oldTime()[cellI]) > SMALL)
		  {
		      t_sol = -gT.oldTime()[cellI]/ddtGT.oldTime()[cellI];
		  }

		  if (t_sol < 0.0)
		  {
		      t_sol = 0.0;
		  }
		  else if (t_sol > runTime.deltaT().value())
		  {
		      t_sol = runTime.deltaT().value();
		  }

                  if(alphaM[cellI]>0.5-SMALL)
		  {  
		  	solidificationTime[cellI] =
		  	runTime.value() - runTime.deltaT().value() + t_sol;
		  	gradTSol[cellI] = gradT[cellI];
		  }
	      }

            if(T[cellI]>Ts.value()) //Reassign values if there is any remelting
	    {
		 solidificationTime[cellI]=-1;
		 gradTSol[cellI] = vector::zero;
	    }
	}
       
	solidificationTime.correctBoundaryConditions();
	gradTSol.correctBoundaryConditions(); 
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	
	if(Power.value() == 0.0)
	{
	    bool solidified = true;
	    if(gMax(alphaL)>1e-6)
	    {
		solidified = false;
	    }

            if (solidified)
            {
                runTime.writeNow();
                Info<< "****All cells solidified****" << endl;
		
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
            	<< nl << endl;
                
		break;
            }
	}
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
