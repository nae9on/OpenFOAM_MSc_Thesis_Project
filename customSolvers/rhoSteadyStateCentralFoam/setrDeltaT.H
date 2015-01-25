/*
{
	// Set the reciprocal time-step from the local Courant number
	surfaceScalarField amaxSfbyDelta
	(
		mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
	);


	rDeltaT.dimensionedInternalField() = max
	(
		1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
		fvc::surfaceSum(mag(amaxSfbyDelta/mesh.magSf()))/(2.0*maxCo)
	);
} 
*/

{    
    // Set the reciprocal time-step from the local Courant number
    Info << "Max. allowed Co: " << maxCo << endl;
    Info << "Max. allowed dT: " << maxDeltaT << endl;

    rDeltaT.dimensionedInternalField() = max
    (   
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(amaxSf)().dimensionedInternalField()
        /(maxCo*mesh.V())
    );
    
    Info << "Local min. dT: " << min(rDeltaT).value() << " : Local max. dT: " << max(rDeltaT).value() << endl;

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();
        
}
