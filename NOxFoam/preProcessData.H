
// Defining dimension for temperature
dimensionedScalar tempDim ("tempDim", dimensionSet(0, 0, 0, 1, 0, 0, 0), 1);

// Taking care of division by zero
dimensionedScalar smallConc ("smallConc", dimensionSet(0, -3, 0, 0, 1, 0, 0), 1e-9);

// Defining gas constant
dimensionedScalar Rgas ("Rgas", dimensionSet(1, 2, -2, -1, -1, 0, 0), 8.3144621);
Info<< "Gas constant " << Rgas << endl << endl;

// Calculating average molecular weight assumning species present are CH4, O2, CO2, H2O, N2
volScalarField M_avg ("M_avg", CH4*M_CH4+O2*M_O2+CO2*M_CO2+H2O*M_H2O+N2*M_N2);

// Calculating average density
volScalarField rhoAvg ("rhoAvg", M_avg*p/(Rgas*T));

// Defining pre-exponential factors
dimensionedScalar ekf1 ("ekf1", dimensionSet(0, 3, -1, 0, -1, 0, 0), 1.8*pow(10,8));
dimensionedScalar ekr1 ("ekr1", dimensionSet(0, 3, -1, 0, -1, 0, 0), 3.8*pow(10,7));
dimensionedScalar ekf2 ("ekf2", dimensionSet(0, 3, -1, -1, -1, 0, 0), 1.8*pow(10,4));
dimensionedScalar ekr2 ("ekr2", dimensionSet(0, 3, -1, -1, -1, 0, 0), 3.81*pow(10,3));
dimensionedScalar ekf3 ("ekf3", dimensionSet(0, 3, -1, 0, -1, 0, 0), 7.1*pow(10,7));
dimensionedScalar ekr3 ("ekr3", dimensionSet(0, 3, -1, 0, -1, 0, 0), 1.7*pow(10,8));
dimensionedScalar eo1  ("eo1", dimensionSet(0, -1.5, 0, 0.5, 0.5, 0, 0), 3.97*pow(10,5));
dimensionedScalar eo2  ("eo2", dimensionSet(0, -1.5, 0, -0.5, 0.5, 0, 0), 36.64);
dimensionedScalar eoh1 ("eoh1", dimensionSet(0, 0, 0, 0.57, 0, 0, 0), 2.129*pow(10,2));

// Defining reaction rates
volScalarField kf1 ("kf1",  ekf1 * exp(-38370*tempDim/T));
volScalarField kr1 ("kr1",  ekr1 * exp(-425*tempDim/T));
volScalarField kf2 ("kf2",  ekf2 * T * exp(-4680*tempDim/T));
volScalarField kr2 ("kr2",  ekr2 * T * exp(-20820*tempDim/T));
volScalarField kf3 ("kf3",  ekf3 * exp(-450*tempDim/T));
volScalarField kr3 ("kr3",  ekr3 * exp(-24560*tempDim/T));
volScalarField ko1 ("ko1",   eo1 * exp(-31090*tempDim/T));
volScalarField ko2 ("ko2",   eo2 * exp(-27123*tempDim/T));
volScalarField koh1("koh1", eoh1 * exp(-4595*tempDim/T));

// Calculating concentration
volScalarField concCH4 ("concCH4", CH4*p/(Rgas*T));
volScalarField concO2 ("concO2",    O2*p/(Rgas*T));
volScalarField concCO2 ("concCO2", CO2*p/(Rgas*T));
volScalarField concH2O ("concH2O", H2O*p/(Rgas*T));
volScalarField concN2 ("concN2",    N2*p/(Rgas*T));
volScalarField concNO ("concNO",    NO*p/(Rgas*T));

// Calculating NO ppm
volScalarField ppmNO ("ppmNO",    NO*pow(10,6));

// Calculating [O], [OH] concentration
////////////////////////////////////////////////////////////////////////////////////////////////////
Info<< "Reading field O if available\n" << endl;
volScalarField O
(
IOobject
(
    "O",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
),

//if (io.headerOk())  { }

// Equilibrium Approach for [O]
//volScalarField concO ("concO",ko1 * pow(T,-0.5) * pow(concO2 , 0.5));
//ko1 * pow(T,-0.5) * pow(concO2 , 0.5) * (Rgas*T)/p

// Partial Equilibrium Approach for [O]
//volScalarField concO ("concO",ko2 * pow(T,0.5) * pow(concO2 , 0.5));
ko2 * pow(T,0.5) * pow(concO2 , 0.5) * (Rgas*T)/p
);
volScalarField concO ("concO", O*p/(Rgas*T));

Info<< "Reading field OH if available\n" << endl;
volScalarField OH
(
IOobject
(
    "OH",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
),
// Partial Equilibrium Approach for [OH]
//volScalarField concOH ("concOH",koh1 * pow(T,-0.57) * pow(concO , 0.5) * pow(concH2O , 0.5));
koh1 * pow(T,-0.57) * pow(concO , 0.5) * pow(concH2O , 0.5) * (Rgas*T)/p
);
volScalarField concOH ("concOH", OH*p/(Rgas*T));

////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate effective diffusivity
volScalarField rhoDTeff ("rhoDTeff", (As * pow(T,1.5) / (T + Ts))/Sc + mut/Sc_t);

// Writitng pre-processed data
concCH4.write();
concO2.write();
concCO2.write();
concH2O.write();
concN2.write();
concNO.write();
concO.write();
concOH.write();

O.write();
OH.write();
M_avg.write();
rhoAvg.write();
rhoDTeff.write();

Info<< "Pre-Processing Completed .................................................................................\n" << endl;
