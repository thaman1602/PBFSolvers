if(sigmaExpMethod == "standard")
{
    sigmaExp = fvc::div
      (
       mu*gradU.T() + lambda*(I*tr(gradD)) - (mu + lambda)*gradD,
       "div(sigma)"
       );
}
else if(sigmaExpMethod == "surface")
{
     sigmaExp = fvc::div
     (
         muf*(mesh.Sf() & fvc::interpolate(gradD.T()))
       + lambdaf*(mesh.Sf() & I*fvc::interpolate(tr(gradD)))
       - (muf + lambdaf)*(mesh.Sf() & fvc::interpolate(gradD))
    );
}
else if(sigmaExpMethod == "decompose")
{
     surfaceTensorField shearGradU = ((I - n*n)&fvc::interpolate(gradD));

     sigmaExp = fvc::div
     (
         mesh.magSf()
        *
         (
            - (muf + lambdaf)*(fvc::snGrad(D)&(I - n*n))
            + lambdaf*tr(shearGradD&(I - n*n))*n
            + muf*(shearGradD&n)
         )
    );
}
else if(sigmaExpMethod == "expLaplacian")
{
    sigmaExp =
       - fvc::laplacian(mu + lambda, solidD, "laplacian(DD,D)")
       + fvc::div
         (
             mu*gradD.T()
           + lambda*(I*tr(gradD)),
             "div(sigma)"
         );
}
else
{
     FatalError << "sigmaExp method " << sigmaExpMethod << " not found!" << endl;
}
