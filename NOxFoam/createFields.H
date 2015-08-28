/*
************************************************* 
Creating variables sourceThermalNO and sourcePromptNO
and assigning them to zero.
*/
    Info<< "Creating volScalarField sourceThermalNO and assigning it to 0\n" << endl;

    volScalarField sourceThermalNO
    (
        IOobject
        (
            "sourceThermalNO",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sourceThermalNO", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    );

    Info<< "Creating volScalarField sourcePromptNO and assigning it to 0\n" << endl;

    volScalarField sourcePromptNO
    (
        IOobject
        (
            "sourcePromptNO",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sourcePromptNO", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    );

//************************************************** Reading p,T,U,mut
    
    Info<< "Reading field p\n" << endl;

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field turbulent viscosity mut\n" << endl;

    volScalarField mut
    (
        IOobject
        (
            "mut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//**************************************************

    Info<< "Reading Properties for NO Computation\n" << endl;

    IOdictionary thermoProperties
    (
        IOobject
        (
            "NO_Properties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading Molar Masses of CH4 O2 CO2 H2O N2 NO from NO_Properties\n" << endl;
    dimensionedScalar M_CH4(thermoProperties.lookup("M_CH4"));
    dimensionedScalar M_O2(thermoProperties.lookup("M_O2"));
    dimensionedScalar M_CO2(thermoProperties.lookup("M_CO2"));
    dimensionedScalar M_H2O(thermoProperties.lookup("M_H2O"));
    dimensionedScalar M_N2(thermoProperties.lookup("M_N2"));
    dimensionedScalar M_NO(thermoProperties.lookup("M_NO"));

    Info<< "Reading Sutherland Transport Coefficients from NO_Properties\n" << endl;
    dimensionedScalar As(thermoProperties.lookup("As"));
    dimensionedScalar Ts(thermoProperties.lookup("Ts"));

    Info<< "Reading Scmdit number and turbulent Scmidt number from NO_Properties\n" << endl;
    dimensionedScalar Sc(thermoProperties.lookup("Sc"));
    dimensionedScalar Sc_t(thermoProperties.lookup("Sc_t"));

    Info<< "Reading flags for thermal and prompt NO from NO_Properties\n" << endl;
    dimensionedScalar flagThermal(thermoProperties.lookup("flagThermal"));
    dimensionedScalar flagPrompt(thermoProperties.lookup("flagPrompt"));

    if (flagThermal.value()) 
    { Info<< "Thermal NOx will be calculated " << endl << endl; }
    else 
    { Info<< "Thermal NOx will not be calculated " << endl << endl; }
    
    if (flagPrompt.value()) 
    { Info<< "Prompt NOx will be calculated " << endl << endl;  }
    else 
    { Info<< "Prompt NOx will not be calculated " << endl << endl; }
    
#   include "createPhi.H"
