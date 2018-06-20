Info << "Reading target field\n" << endl;
volScalarField Td
(
    IOobject
    (
	"Td",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);
