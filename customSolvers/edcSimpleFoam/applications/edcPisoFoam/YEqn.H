{
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];


    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix Yeqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                chemistry.RR(i) //added
            );
            
            eqnResidual = Yeqn.solve(mesh.solver("Yi")).initialResidual();
            maxResidual = max(eqnResidual, maxResidual);            



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
