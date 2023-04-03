/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description
    TBD

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "wallDist.H"
#include "fvCFD.H"
#include "fvMeshFunctionObject.H"
#include "searchableSurfaces.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "triSurfaceSearch.H"
#include "IOobjectList.H"


using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    const word dictName("labelRegionDict");

    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary dict(dictIO);

    Info<< "Time now = " << runTime.timeName() << nl <<endl;

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const dictionary& labelsDict = dict.subDict("labels");
    word labelsOutputField = word(labelsDict.lookup("outputField"));
    word labelsSourceField = word(labelsDict.lookup("sourceField"));
    scalar labelsInletID = readScalar(labelsDict.lookup("inletID"));
    scalar labelsOutletID = readScalar(labelsDict.lookup("outletID"));
    scalar labelsWallID = readScalar(labelsDict.lookup("wallID"));

    Info<< "Writing "<< labelsOutputField << " field" << nl <<endl;

    wordRes labelsInletPatchNames =
        labelsDict.lookupOrDefault<wordRes>("inletPatches", wordRes());
    wordRes labelsOutletPatchNames =
        labelsDict.lookupOrDefault<wordRes>("outletPatches", wordRes());
    wordRes labelsWallPatchNames =
        labelsDict.lookupOrDefault<wordRes>("wallPatches", wordRes());

    volScalarField sourceField
    (
        IOobject
        (
            labelsSourceField,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    sourceField.rename(labelsOutputField);

    labelList labelsPatchInletIDs
    (
        pbm.patchSet(labelsInletPatchNames).sortedToc()
    );

    labelList labelsPatchOutletIDs
    (
        pbm.patchSet(labelsOutletPatchNames).sortedToc()
    );

    labelList labelsPatchWallIDs
    (
        pbm.patchSet(labelsWallPatchNames).sortedToc()
    );

    const fvPatchList& patches = mesh.boundary();

    forAll(labelsPatchInletIDs, labelsPatchInletID) {
        const fvPatch& currPatch =
            patches[labelsPatchInletIDs[labelsPatchInletID]];
        forAll(currPatch, facei)
        {
             label faceCelli = currPatch.faceCells()[facei];
             sourceField[faceCelli] = labelsInletID;
        }
    }

    forAll(labelsPatchOutletIDs, labelsPatchOutletID) {
        const fvPatch& currPatch =
            patches[labelsPatchOutletIDs[labelsPatchOutletID]];
        forAll(currPatch, facei)
        {
             label faceCelli = currPatch.faceCells()[facei];
             sourceField[faceCelli] = labelsOutletID;
        }
    }

    forAll(labelsPatchWallIDs, labelsPatchWallID) {
        const fvPatch& currPatch =
            patches[labelsPatchWallIDs[labelsPatchWallID]];
        forAll(currPatch, facei)
        {
             label faceCelli = currPatch.faceCells()[facei];
             sourceField[faceCelli] = labelsWallID;
        }
    }

    sourceField.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
