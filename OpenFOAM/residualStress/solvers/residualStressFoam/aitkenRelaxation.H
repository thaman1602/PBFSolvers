// aitken acceleration

// update delta
aitkenDelta.storePrevIter();
aitkenDelta = (D - D.prevIter()) / aitkenInitialRes;

// update relaxation factor
if(iCorr == 0)
{
    aitkenTheta = 0.01;
//     if(mesh.relax(D.name()))
//       {
//     aitkenTheta = mesh.relaxationFactor(D.name());
//       }
}
else
{
    vectorField b = aitkenDelta.internalField() - aitkenDelta.prevIter().internalField();
    // scalar sumMagB = gSum(mag(b));
    scalar sumMagB = gSum(magSqr(b));
    if(sumMagB < SMALL)
    {
        // Warning << "Aitken under-relaxation: denominator less then SMALL"
        //    << endl;
        sumMagB += SMALL;
    }

    aitkenTheta = -aitkenTheta*
        gSum(aitkenDelta.prevIter().internalField() & b)/sumMagB;
}

// correction to the latest D
D += aitkenTheta*aitkenDelta*aitkenInitialRes;
