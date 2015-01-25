Info<< nl << "Reading thermophysical properties\n" << endl;

autoPtr<edcChemistryModel> edcChemistry
(
    edcChemistryModel::New(mesh)
);

edcChemistryModel& chemistry = edcChemistry();
hsReactionThermo& thermo = chemistry.thermo();

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

//access reaction data; re-cast thermo as reactionMixture
const PtrList<Reaction<gasThermoPhysics> >& reactions =
            dynamic_cast<const reactingMixture<gasThermoPhysics>&>  (thermo);

//access specie thermo data 
const PtrList<gasThermoPhysics> & specieThermo =
            dynamic_cast<const reactingMixture<gasThermoPhysics>&>  (thermo).speciesData();    


word inertSpecie(thermo.lookup("inertSpecie"));    

Switch limitTemperature(thermo.lookup("limitTemperature"));
dimensionedScalar TMin(thermo.lookup("TMin"));
dimensionedScalar TMax(thermo.lookup("TMax"));    

bool additionalOutput=readBool(chemistry.lookup("additionalOutput"));


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

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
    

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll (Y, i)
{
    fields.add(Y[i]);
}
//fields.add(hs);

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
