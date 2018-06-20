Info << "Reading physicalProperties\n" << endl;
IOdictionary physicalProperties
(
    IOobject
    (
	"physicalProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

dimensionedScalar k(physicalProperties.lookup("k"));
dimensionedScalar beta(physicalProperties.lookup("beta"));

Info << "Reading numericalParameters\n" << endl;
IOdictionary numericalParameters
(
    IOobject
    (
	"numericalParameters",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

scalar gamma(readScalar(numericalParameters.lookup("gamma")));
label incount(readLabel(numericalParameters.lookup("incount")));
wordList movingPatches(numericalParameters.lookup("movingPatches"));

Info << "Reading field T\n" << endl;
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

surfaceScalarField snGradT
(
    IOobject
    (
	"snGradT",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    fvc::snGrad(T)
);

Info << "Calculating heat flux density q\n" << endl;
volVectorField q
(
    IOobject
    (
	"q",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    -k*(fvc::grad(T))
);

Info << "Reading field f\n" << endl;
volScalarField f
(
    IOobject
    (
	"f",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);

#include "readTarget.H"

/*
Info << "Reading target field\n" << endl;
volScalarField Td
(
    IOobject
    (
	"Td",
	runTime.timeName(),
	mesh,
	IOobject::MUST_READ,
	IOobject::AUTO_WRITE
    ),
    mesh
);
*/

Info << "Reading adjoint field\n" << endl;
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

surfaceScalarField snGradp
(
    IOobject
    (
	"snGradp",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    fvc::snGrad(p)
);

surfaceScalarField sensitivity
(
    IOobject
    (
	"sensitivity",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    fvc::interpolate( 0.5 * (T - Td) * (T - Td) ) + k*snGradT*snGradp
);

scalarField volField = mesh.V();

