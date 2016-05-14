    const dictionary& fluidPimple = mesh.solutionDict().subDict("PIMPLE");

    const int nCorr =
        fluidPimple.lookupOrDefault<int>("nCorrectors", 1);

    const int nNonOrthCorr =
        fluidPimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        fluidPimple.lookupOrDefault("momentumPredictor", true);
