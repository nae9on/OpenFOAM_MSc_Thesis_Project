/* ------------------- furnaceFoam source code */

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiCombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "radiationModel.H"

int main(int argc, char *argv[])
{
    // 1. Pre-processing
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
        // 2. Calculate T
        {
         // Read the control parameters used by setDeltaT
         #include "readTimeControls.H" 

         // Calculates and outputs the mean and maximum Courant Numbers.
         #include "compressibleCourantNo.H" 

         /* Reset the timestep to maintain a constant maximum courant Number.
         Reduction of time-step is immediate, but increase is damped to avoid
         unstable oscillations. */
         #include "setDeltaT.H" 
        
         runTime++;
         Info<< "Time = " << runTime.timeName() << nl << endl;
        }
        
        // 3. Solve the continuity equation to compute density.
        #include "rhoEqn.H" 

        while (pimple.loop()) // outer iteration nOuterCorrectors = 1 (= 1 means PISO mode)
        {
            // 4. Momentum Predictor
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
            // if it is on once the solve method is executed U is updated using the old
            //  pressure field.
	    if (pimple.momentumPredictor()) // set No in PIMPLE control
	    {
		solve(UEqn == -fvc::grad(p));

		fvOptions.correct(U);
		K = 0.5*magSqr(U);
	    }
            // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


            // 5. Species Transport
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
            
            // 6. Energy Transport
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

                    // 7. Radiation Transport
                    // Solve transport equation for G
		    radiation->correct(); 
                    
                    // 8. Update thermophysical properties
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
