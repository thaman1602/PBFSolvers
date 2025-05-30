t = mesh.time();

Power.value() = laserPowerSeries(runTime.time().timeOutputValue());
laserPos.value() = laserPositionSeries(runTime.time().timeOutputValue());

dimensionedScalar x = laserPos.component(0);
dimensionedScalar y = laserPos.component(1);

volScalarField gaussian = Foam::exp(-2.0*(pow(x - px,2)+pow(y - py,2))/(r0*r0));

#   include "updateFLB.H"

laserFlux = 2.0 * flb * kappaM * Power * gaussian/(M_PI*r0*r0);
