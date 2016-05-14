    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& mainPimpleDict = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        mainPimpleDict.lookupOrDefault<int>("nOuterCorrectors", 1);
