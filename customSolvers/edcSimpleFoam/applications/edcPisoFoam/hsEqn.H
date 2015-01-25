{
    fvScalarMatrix hsEqn
    (
       fvm::ddt(rho, hs)
      +  fvm::div(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)

      -  radiation->Shs(thermo)  //heat source sink due to radiation
     -  chemistrySh     // heat release from reaction
     ==
        DpDt
    );

    hsEqn.relax();

    eqnResidual = hsEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);


    if (limitTemperature)
    {
        //Limit temperature to TMin and TMax

        // define temporary fields for hMin and hMax as copy of h
        tmp<volScalarField> hsMin
            (new volScalarField
               (
                 IOobject
                 (
                     "hsMin",
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 hs
               )
           );
        tmp<volScalarField> hsMax
            (new volScalarField
               (
                 IOobject
                 (
                     "hsMax",
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE
                 ),
                 hs
               )
           );


        gasThermoPhysics mixture(specieThermo[0]);


        forAll(hs.internalField(), cellI)
        {
            mixture = Y[0].internalField()[cellI]/specieThermo[0].W()*specieThermo[0];

            for (label specieI=1; specieI < composition.Y().size(); specieI++)        
            {
        	    mixture += Y[specieI].internalField()[cellI]/specieThermo[specieI].W()*specieThermo[specieI];
            }    
            hsMin().internalField()[cellI]=mixture.Hs(TMin.value());
            hsMax().internalField()[cellI]=mixture.Hs(TMax.value());

        }

        forAll(hs.boundaryField(), patchI)
        {
            
            forAll(hs.boundaryField()[patchI], faceI)
            {
                mixture= Y[0].boundaryField()[patchI][faceI]/specieThermo[0].W()*specieThermo[0];
                for (label specieI=1; specieI < composition.Y().size(); specieI++)        
                {
                    mixture+= Y[specieI].boundaryField()[patchI][faceI]/specieThermo[specieI].W()*specieThermo[specieI];
                }
                hsMin().boundaryField()[patchI][faceI]=mixture.Hs(TMin.value());
                hsMax().boundaryField()[patchI][faceI]=mixture.Hs(TMax.value());

            }
        }

      



        Info << "Limiting temperature between " << TMin << " and " << TMax << endl;
        //bound h to TMin
        hs=max(hs, hsMin());
        //bound h to TMax
        hs=min(hs, hsMax());

    }

    thermo.correct();
    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;

    radiation->correct();
    
}
