/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionFoam

Description
    Combination of heatConductionFoam and buoyantFoam for conjugate heat
    transfer between solid regions and fluid regions. Both regions include
    the fvOptions framework.

    It handles secondary fluid or solid circuits which can be coupled
    thermally with the main fluid region. i.e radiators, etc.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvIOoptionList.H"
#include "coordinateSystem.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    #include "createFluidMeshes.H"
    #include "createSolidMeshes.H"
    
    pimpleControl pimple(fluidRegions[0]);

    #include "createFluidFields.H"
    #include "createSolidFields.H"

    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"


    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    while (runTime.run())
    {
        #include "createTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }


        // --- PIMPLE loop
//        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        while (pimple.loop())
        {
//            bool finalIter = oCorr == nOuterCorr-1;

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
//                #include "solveFluid.H"
//                if (finalIter)
//                {
//                    mesh.data::add("finalIteration", true);
//                }

                if (frozenFlow)
                {
                    #include "EEqn.H"
                }
                else
                {
                    if (pimple.firstIter())
                    {
                        #include "rhoEqn.H"
                    }

                    #include "UEqn.H"
                    #include "EEqn.H"

                    // --- PISO loop
//                    for (int corr=0; corr<nCorr; corr++)
                    while (pimple.correct())
                    {
                        #include "pEqn.H"
                    }
                    if (pimple.turbCorr())
                    {
                        turb.correct();
                    }

                    rho = thermo.rho();
                }
    
//                if (finalIter)
//                {
//                    mesh.data::remove("finalIteration");
//                }
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
//                #include "solveSolid.H"
//                if (finalIter)
//                {
//                    mesh.data::add("finalIteration", true);
//                }

//                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                while (pimple.correctNonOrthogonal())
                {
                    tmp<fvScalarMatrix> hEqn
                    (
                        fvm::ddt(betav*rho, h)
                      - (
                           thermo.isotropic()
                         ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
                         : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
                        )
                      ==
                        fvOptions(rho, h)
                    );

                    hEqn().relax();

                    fvOptions.constrain(hEqn());

            //        hEqn().solve(mesh.solver(h.select(finalIter)));
                    hEqn().solve(mesh.solver(h.select(pimple.finalIter())));

                    fvOptions.correct(h);
                }

                thermo.correct();

                Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                    << max(thermo.T()).value() << endl;

//                if (finalIter)
//                {
//                    mesh.data::remove("finalIteration");
//                }
            }

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
