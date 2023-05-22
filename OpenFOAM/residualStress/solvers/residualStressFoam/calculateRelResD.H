{
    // displacement
    scalar maxDD = gMax(mag(D.internalField() - D.oldTime().internalField()));
    relResD =
        gMax
        (
            mag(D.internalField() - D.prevIter().internalField())
           /(maxDD + SMALL)
        );
}
