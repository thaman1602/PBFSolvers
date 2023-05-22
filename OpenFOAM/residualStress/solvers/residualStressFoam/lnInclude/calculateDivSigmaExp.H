if(divSigmaExpMethod == "standard")
{
    divSigmaExp = fvc::div
    (
        mu*gradD.T() + lambda*(I*tr(gradD)) - (mu + lambda)*gradD,
        "div(sigma)"
    )
    - gradThreeKalphaDeltaT;
    Info<<"min/max(divSigExp)"<<max(divSigmaExp).value()<<", "<<min(divSigmaExp).value()<<endl;
}
else if(divSigmaExpMethod == "surface")
{
    divSigmaExp = fvc::div
    (
        muf*(mesh.Sf() & fvc::interpolate(gradD.T()))
      + lambdaf*(mesh.Sf() & I*fvc::interpolate(tr(gradD)))
      - (muf + lambdaf)*(mesh.Sf() & fvc::interpolate(gradD))
      - threeKalphaDeltaTf
    );
}
else if(divSigmaExpMethod == "decompose")
{
    snGradD = fvc::snGrad(D);
    surfaceTensorField shearGradD
    (
        IOobject
        (
            "shearGrad(D)",
            runTime.timeName(),
            solidMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        solidMesh,
        dimensionedTensor("zero", dimless, tensor::zero)
     );

    //surfaceTensorField shearGradD = ((I - n*n) & fvc::interpolate(gradD));
    shearGradD = ((I - n*n) & fvc::interpolate(gradD));

    divSigmaExp = fvc::div
    (
        solidMesh.magSf()*
        (
          - (muf + lambdaf)*(snGradD & (I - n*n))
          + lambdaf*tr(shearGradD & (I - n*n))*n
          + muf*(shearGradD & n)
        )
      - threeKalphaDeltaTf
    );
}
/*
else if(divSigmaExpMethod == "expLaplacian")
{
    divSigmaExp =
      - fvc::laplacian(mu + lambda, U, "laplacian(DU,U)")
      + fvc::div(mu*gradU.T() + lambda*(I*tr(gradU)), "div(sigma)");
}
*/
else
{
    FatalErrorIn(args.executable())
        << "divSigmaExp method " << divSigmaExpMethod << " not found!" << endl;
}
