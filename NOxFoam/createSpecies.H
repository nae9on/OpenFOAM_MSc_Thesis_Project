
//************************************************** Reading Species

    Info<< "Reading field CH4\n" << endl;

    volScalarField CH4
    (
        IOobject
        (
            "CH4",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field O2\n" << endl;

    volScalarField O2
    (
        IOobject
        (
            "O2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field CO2\n" << endl;

    volScalarField CO2
    (
        IOobject
        (
            "CO2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field H2O\n" << endl;

    volScalarField H2O
    (
        IOobject
        (
            "H2O",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field N2\n" << endl;

    volScalarField N2
    (
        IOobject
        (
            "N2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field NO\n" << endl;

    volScalarField NO
    (
        IOobject
        (
            "NO",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //NO = NO/2;
    //NO.write();
