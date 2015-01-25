    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<edcChemistryModel> edcChemistry
    (
        edcChemistryModel::New(mesh)
    );
    edcChemistryModel& chemistry = edcChemistry();



    hsReactionThermo& thermo = chemistry.thermo();

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();
    

    //access specie thermo data 
    const PtrList<gasThermoPhysics> & specieThermo =
             dynamic_cast<const reactingMixture<gasThermoPhysics>&>  (thermo).speciesData();    
    

    word inertSpecie(thermo.lookup("inertSpecie"));    
    Switch limitTemperature(thermo.lookup("limitTemperature"));
    dimensionedScalar TMin(thermo.lookup("TMin"));
    dimensionedScalar TMax(thermo.lookup("TMax"));    
    

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    
    


    volScalarField& p = thermo.p();
    volScalarField& hs = thermo.hs();
    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();
    
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

#   include "compressibleCreatePhi.H"

    surfaceScalarField UBlendingFactor
    (
        IOobject
        (
            "UBlendingFactor",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(Y[0])

    );    
    UBlendingFactor = 0.5;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    dimensionedScalar pMin
    (
        mesh.solutionDict().subDict("SIMPLE").lookup("pMin")
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );
    dimensionedScalar initialMass = fvc::domainIntegrate(rho);
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll (Y, i)
    {
        fields.add(Y[i]);
    }
//    fields.add(hs);    
    
DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.)
);    
