#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "radiationModel.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"

    #include "createFvOptions.H"

    #include "createRadiationModel.H"

    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    
    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run()) // to iterate in time.
    {
        #include "readTimeControls.H" // Read the control parameters used by setDeltaT
        #include "compressibleCourantNo.H" // Calculates and outputs the mean and maximum Courant Numbers.

        /* Reset the timestep to maintain a constant maximum courant Number.
        Reduction of time-step is immediate, but increase is damped to avoid
        unstable oscillations. */
        #include "setDeltaT.H" 


        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "rhoEqn.H" // Solve the continuity equation for density.

        while (pimple.loop()) // outer iteration nOuterCorrectors = 1 (= 1 means PISO mode)
        {
            // #include "UEqn.H"
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	    fvVectorMatrix UEqn
	    (
		fvm::ddt(rho, U)
	      + fvm::div(phi, U)
	      + turbulence->divDevRhoReff(U)
	     ==
		rho*g
	      + fvOptions(rho, U)
	    );

	    UEqn.relax();

	    fvOptions.constrain(UEqn);

            // if momentum predictor is turned off, U from previous time-step is used, 
            // if it is on once the solve method is executed U is updated using the old pressure field.
	    if (pimple.momentumPredictor()) // set No in PIMPLE control
	    {
		solve(UEqn == -fvc::grad(p));

		fvOptions.correct(U);
		K = 0.5*magSqr(U);
	    }
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

            // #include "YEqn.H"
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
		tmp<fv::convectionScheme<scalar> > mvConvection
		(
		    fv::convectionScheme<scalar>::New
		    (
			mesh,
			fields,
			phi,
			mesh.divScheme("div(phi,Yi_h)")
		    )
		);

                // Solve transport equation for individual species mass fractions
		{
		    reaction->correct();
		    dQ = reaction->dQ();
		    label inertIndex = -1;
		    volScalarField Yt(0.0*Y[0]);

		    forAll(Y, i)
		    {
			if (Y[i].name() != inertSpecie)
			{
			    volScalarField& Yi = Y[i];

			    fvScalarMatrix YiEqn
			    (
				fvm::ddt(rho, Yi)
			      + mvConvection->fvmDiv(phi, Yi)
			      - fvm::laplacian(turbulence->muEff(), Yi)
			     ==
				reaction->R(Yi)
			      + fvOptions(rho, Yi)
			    );

			    YiEqn.relax();

			    fvOptions.constrain(YiEqn);

			    YiEqn.solve(mesh.solver("Yi"));

			    fvOptions.correct(Yi);

			    Yi.max(0.0);
			    Yt += Yi;
			}
			else
			{
			    inertIndex = i;
			}
		    }

		    Y[inertIndex] = scalar(1) - Yt;
		    Y[inertIndex].max(0.0);
		}
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

            // #include "EEqn.H"
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
                // Solve energy equation to obtain new temperature field
		{
		    volScalarField& he = thermo.he();

		    fvScalarMatrix EEqn
		    (
			fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
		      + fvc::ddt(rho, K) + fvc::div(phi, K)
		      + (
			    he.name() == "e"
			  ? fvc::div
			    (
				fvc::absolute(phi/fvc::interpolate(rho), U),
				p,
				"div(phiv,p)"
			    )
			  : -dpdt
			)
		      - fvm::laplacian(turbulence->alphaEff(), he)
		     ==
			reaction->Sh()        // Source term for combustion
		      + radiation->Sh(thermo) // Source term for radiation
		      + fvOptions(rho, he)
		    );

		    EEqn.relax();

		    fvOptions.constrain(EEqn);

		    EEqn.solve();

		    fvOptions.correct(he);

		    radiation->correct(); // Solve transport equation for G

		    thermo.correct();

		    Info<< "min/max(T) = "
			<< min(T).value() << ", " << max(T).value() << endl;
		}
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

            // --- Pressure corrector loop
            while (pimple.correct()) // - inner iteration nCorrectors = 2
            {
                #include "pEqn.H"
            }

            // Its function is to determine whether to correct turbulence at every outer 
            // iteration or only at the last iteration.
            // turbOnFinalIterOnly is set by default to true
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        } // Done PIMPLE - outer iteration

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    } // Done runTime.run()

    Info<< "End\n" << endl;

    return 0;
} // Done - main
