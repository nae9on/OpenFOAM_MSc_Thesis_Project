{

    /* rewritten */
    chemistry.solve(0.,10.);
    if(runTime.write())
    {
        forAll(Y,specieI)
        {
            chemistry.RR(specieI)->write();
        }
    }

    chemistrySh = chemistry.Sh();
  
}
