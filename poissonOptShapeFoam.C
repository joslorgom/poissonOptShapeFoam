/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    poissonOptShapeFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "dynamicFvMesh.H"
#include "primitivePatchInterpolation.H"
#include "fixedValuePointPatchField.H"

#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"

    // Disable primal and adjoint solvers output
    lduMatrix::debug = 0;
    solverPerformance::debug = 0;

    scalar volSum = gSum(volField);
    Info << "Initial Volume is " << volSum << nl << endl;

    scalar Ja = VGREAT;
    scalar J = VGREAT;
    scalar dJ = 0;
    label count = 0;

    std::ofstream file("results.csv");
    file.precision(16);

    while (runTime.loop() && count < incount)
    {

	// Primal equation
	solve(fvm::laplacian(k, T) + f);
	// Heat flux density
	q = -k*(fvc::grad(T));
	// Normal derivative of temperature field at all surfaces
	snGradT = fvc::snGrad(T);

	#include "readTarget.H"

	// Adjoint equation
	solve(fvm::laplacian(k, p) + T - Td);
	// Normal derivative of adjoint field at all surfaces
	snGradp = fvc::snGrad(p);

	volField = mesh.V();
	volSum = gSum(volField);

	// Update cost function value
	Ja = J;
	J = 0.5 * gSum( volField * (T.internalField() - Td.internalField()) * (T.internalField() - Td.internalField()) );
	dJ = J - Ja;

	if(dJ > 0)
	    {count ++;}
	else
	    {count = 0;}

	Info << setprecision(8);
	Info << nl << "Iteration " << runTime.timeName() << " - Cost " << J << " - Mesh Volume " << volSum << nl << endl;

	file << runTime.value() << "," << J << "," << volSum << nl;

	// Sensitivity field
	sensitivity = fvc::interpolate( 0.5*(T - Td)*(T - Td) ) + k*snGradT*snGradp;

	// Move mesh
	forAll(movingPatches,i)
	{
	    const word& patchName = movingPatches[i];
	    label patchID = mesh.boundaryMesh().findPatchID(patchName);

	    const polyPatch& cPatch = mesh.boundaryMesh()[patchID]; 

	    scalar thetaSum(0);
	    scalar areaSum(0);

	    // Enforce constant volume
	    forAll(cPatch,faceI) 
	    {
		thetaSum = thetaSum + sensitivity.boundaryField()[patchID][faceI] * mesh.magSf().boundaryField()[patchID][faceI];
		areaSum = areaSum + mesh.magSf().boundaryField()[patchID][faceI];
	    }
	    sensitivity = sensitivity - dimensionedScalar("r", sensitivity.dimensions(), ( thetaSum / areaSum ) );

	    pointVectorField& PointDisplacement = const_cast<pointVectorField&>
	    (
	    	mesh.objectRegistry::lookupObject<pointVectorField>
	    	(
		    "pointDisplacement"
	  	)
	    );

	    //Get the vector field of the patch
	    vectorField &pDisp = refCast<vectorField>(PointDisplacement.boundaryFieldRef()[patchID]);

	    //Find the relevant size of the vector and declare a vectorField.
	    int Psize= pDisp.size();
	    vectorField dispVals(Psize);

	    // displacement based on the external calculation
	    scalarField sensitivityPatch = sensitivity.boundaryField()[patchID]/gMax(mag(sensitivity.boundaryField()[patchID]));

	    //- set-up interpolator
	    primitivePatchInterpolation patchInterpolator( mesh.boundaryMesh()[patchID] );

	    //- perform interpolation
	    scalarField faceValues = patchInterpolator.faceToPointInterpolate(sensitivityPatch);

	    vectorField &PointPointer = refCast<vectorField>(PointDisplacement.boundaryFieldRef()[patchID]);
	    vectorField PointNormalVector = mesh.boundaryMesh()[patchID].pointNormals();

	    // loop over points to move the nodes
	    forAll(dispVals, index)
	    {
	    	dispVals[index].x() = PointPointer[index].x() - gamma * faceValues[index] * PointNormalVector[index].x();
            	dispVals[index].y() = PointPointer[index].y() - gamma * faceValues[index] * PointNormalVector[index].y();
            	dispVals[index].z() = PointPointer[index].z() - gamma * faceValues[index] * PointNormalVector[index].z();
	    }

	    PointDisplacement.boundaryFieldRef()[patchID] == dispVals;
	}

	mesh.update();

	runTime.write();
    }

    file.close();

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
