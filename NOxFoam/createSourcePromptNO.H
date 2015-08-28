
//scalar n_p = 1;
//scalar phi_p = 1;

Info<< "Reading parameters for prompt NOx from NO_Properties\n" << endl;
dimensionedScalar n_p(thermoProperties.lookup("n_p"));
dimensionedScalar phi_p(thermoProperties.lookup("phi_p")); 

dimensionedScalar f_p ("f_p", 4.75 + 0.0819*n_p - 23.2*phi_p + 32*pow(phi_p,2) -12.2*pow(phi_p,3));

// Calculating oxygen reaction order
volScalarField o2MoleFrac ("o2MoleFrac", O2*M_avg/M_O2);

volScalarField a_p
(
IOobject
(
    "a_p",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
),
mesh,
scalar(1)
);


volScalarField o2_a_p
(
IOobject
(
    "o2_a_p",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
),
mesh,
scalar(1)
);


volScalarField kpr_temp
(
IOobject
(
    "kpr_temp",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
),
mesh,
scalar(1)
);

forAll (o2MoleFrac, celli)
{
    if (o2MoleFrac[celli] <= 4.1e-03) {a_p[celli] = 1;}
    else if ( o2MoleFrac[celli] > 4.1e-03 && o2MoleFrac[celli] <= 1.11e-02) {a_p[celli] = -3.95-0.9*Foam::log (o2MoleFrac[celli]);}
    else if ( o2MoleFrac[celli] > 1.11e-02 && o2MoleFrac[celli] < 0.03) {a_p[celli] = -0.35-0.1*Foam::log (o2MoleFrac[celli]);}
    else {a_p[celli] = 0;}

    o2_a_p[celli] = Foam::pow(concO2[celli],a_p[celli]);
    kpr_temp[celli] = Foam::pow(Rgas.value()*T[celli]/p[celli],a_p[celli]+1);
}

//All dimensions will be adjusted here
dimensionedScalar ekpr ("ekpr", dimensionSet(0, 3, -1, 0, -1, 0, 0), 6.4*pow(10,6));

volScalarField kpr ("kpr", ekpr*kpr_temp);

dimensionedScalar Ea_p ("Ea_p", dimensionSet(1, 2, -2, 0, -1, 0, 0), 303474.125);

sourcePromptNO = M_NO * f_p * kpr * o2_a_p * concN2 * concCH4 * exp(-Ea_p/(Rgas*T));

Info<< "Source Prompt NO Created .................................................................................\n" << endl;

