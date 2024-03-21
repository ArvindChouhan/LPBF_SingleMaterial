#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar cps("cps", dimensionSet(0, 2, -2, -1, 0, 0, 0), transportProperties);
dimensionedScalar cpl("cpl", dimensionSet(0, 2, -2, -1, 0, 0, 0), transportProperties);
dimensionedScalar cp2("cp2", dimensionSet(0, 2, -2, -1, 0, 0, 0), transportProperties);

dimensionedScalar ks("ks", dimensionSet(1, 1, -3, -1, 0, 0, 0), transportProperties);
dimensionedScalar kl("kl", dimensionSet(1, 1, -3, -1, 0, 0, 0), transportProperties);
dimensionedScalar k2("k2", dimensionSet(1, 1, -3, -1, 0, 0, 0), transportProperties);

dimensionedScalar beta1("beta1", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);
dimensionedScalar beta2("beta2", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);

dimensionedScalar Ts("Ts", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar Tl("Tl", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar Tref("Tref", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar L("L", dimensionSet(0, 2, -2, 0, 0, 0, 0), transportProperties);
dimensionedScalar C("C", dimensionSet(1, -3, -1, 0, 0, 0, 0), transportProperties);
dimensionedScalar b("b", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar dsigmadt("dsigmadt", dimensionSet(1, 0, -2, -1, 0, 0, 0), transportProperties);
dimensionedScalar E("E", dimensionSet(0, 0, 0, 0, 0, 0, 0), transportProperties);
dimensionedScalar h1("h1", dimensionSet(1, 0, -3, -1, 0, 0, 0), transportProperties);
dimensionedScalar sig("sig", dimensionSet(1, 0, -3, -4, 0, 0, 0), transportProperties);
dimensionedScalar P_0("P_0", dimensionSet(1, -1, -2, 0, 0, 0, 0), transportProperties);
dimensionedScalar L_v("L_v", dimensionSet(0, 2, -2, 0, 0, 0, 0), transportProperties);
dimensionedScalar M("M", dimensionSet(1, 0, 0, 0, -1, 0, 0), transportProperties);
dimensionedScalar T_v("T_v", dimensionSet(0, 0, 0, 1, 0, 0, 0), transportProperties);
dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), transportProperties);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField lf
(
    IOobject
    (
        "lf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
     max(min((T-Ts)/(Tl-Ts),scalar(1)),scalar(0)) 
);

volScalarField k = alpha1*(ks*(scalar(1)-lf)+kl*lf) +  alpha2*k2;


volScalarField rhoCp = alpha1*rho1*(cps*(scalar(1)-lf)+cpl*lf) +  alpha2*rho2*cp2;
surfaceScalarField rhoCpphi = fvc::interpolate(rhoCp)*phi;

volScalarField beta = alpha1*beta1 +  alpha2*beta2;
volScalarField A = -C * alpha1 * sqr(scalar(1)-lf)/ (pow(lf,scalar(3))+b) ;

// latent heat source terms (implicit)
volScalarField lfPRELOOP = lf;
volScalarField lfold = lf;
dimensionedScalar dFdT = 1 / (Tl - Ts);
volScalarField invF = lfold * (Tl - Ts) + Ts;
volScalarField dH = alpha1*rho * L * pos(Tl - T) * pos(T - Ts); //alpha1 multiplication
volScalarField Sp = - dH/runTime.deltaT() * dFdT;
volScalarField Sc = dH/runTime.deltaT() * (lfPRELOOP - lf + dFdT * invF);

volScalarField kappa = mixture.Arv();

///////////////////////////////////////////////////////////////////////////////// Radiation Model /////////////////////////////////////////////////////////////////////////////////////////////


volScalarField alphaT = alpha1;

autoPtr<radiation::radiationModel> radiation
(
        radiation::radiationModel::New(alphaT)
);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






