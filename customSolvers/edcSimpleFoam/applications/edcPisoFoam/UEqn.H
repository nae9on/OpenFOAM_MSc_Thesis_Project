    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     == 
        rho*g
    );
    eqnResidual = UEqn.solve().initialResidual();
    UEqn.relax();
    maxResidual = max(eqnResidual, maxResidual);
    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

