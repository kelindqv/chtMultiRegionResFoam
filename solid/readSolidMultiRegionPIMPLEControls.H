const dictionary& solidPimple = mesh.solutionDict().subDict("PIMPLE");

int nNonOrthCorr =
    solidPimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
