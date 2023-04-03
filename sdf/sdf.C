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


autoPtr<triSurface> loadSurface
(
    const Foam::Time& runTime,
    const fileName& surfName,
    const scalar scaleFactor
)
{
    Info<< "Reading surface " << surfName << nl;
    if (scaleFactor > 0)
    {
        Info<<"Scaling : " << scaleFactor << nl;
    }

    const fileName fallback =
        runTime.constantPath()/triSurfaceMesh::meshSubDir/surfName;

    autoPtr<triSurface> surfPtr;
    if (isFile(surfName))
    {
        surfPtr.set(new triSurface(surfName, scaleFactor));
    }
    else if (isFile(fallback))
    {
        surfPtr.set(new triSurface(fallback, scaleFactor));
    }
    else
    {
        FatalErrorInFunction
            << "No such file:" << surfName << exit(FatalError);
    }

    return surfPtr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    /*argList::addOption
    (
        "patches",
        "wordRes",
        "Specify a patch or list of patches"
        "Eg, 'bottomWall' or '(bottomWall topWall \"wall.*\")'"
    );

    argList::addOption
    (
        "surface",
        "name",
        "Specify a surface name for a file in constant/triSurface"
        "Eg, wall.stl"
    );*/

    argList::addNote
    (
        "Calculates signed distance fields"
    );
 
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    Info<< timeDirs <<endl;

    const word dictName("sdfDict");

    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    IOdictionary dict(dictIO);

    Info<< "Time now = " << runTime.timeName() << nl <<endl;

    // 1 - Surface
    const dictionary& surfaceDict = dict.subDict("surface");
    const dictionary& geometryDict = surfaceDict.subDict("geometry");
    word surfaceOutputField = word(surfaceDict.lookup("outputField"));

    Info<< "Writing "<< surfaceOutputField << " field" << nl <<endl;

    wordList surfNames = geometryDict.toc();
    fileName surfName = surfNames[0];
    triSurface surf = loadSurface(runTime, surfName, -1)();
    triSurfaceSearch query(surf);

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    boolList isInside = query.calcInside(mesh.cellCentres());

    const searchableSurfaces searchSurf
    (
        IOobject
        (
            "abc", // dummy
            mesh.time().constant(),
            "triSurface",
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        true
    );

    volScalarField distance
    (
        IOobject
        (
             surfaceOutputField,
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, Zero)
    );

    const pointField& cc = mesh.C();

    labelList surfaces;
    List<pointIndexHit> nearestInfo;
    searchSurf.findNearest
    (
         cc,
         scalarField(cc.size(), GREAT),
         surfaces,
         nearestInfo
    );

    forAll(nearestInfo, celli)
    {
         scalar magDistance = mag(nearestInfo[celli].hitPoint()-cc[celli]);
         if (isInside[celli] == 1)
             distance[celli] = -magDistance;
         else
             distance[celli] = magDistance;
    }
    distance.correctBoundaryConditions();

    distance.write();

    // 2 - Patches
    const dictionary& patchesDict = dict.subDict("patches");
    word patchesOutputField = word(patchesDict.lookup("outputField"));

    Info<< "Writing "<< patchesOutputField << " field" << nl <<endl;

    wordRes patchNames =
        patchesDict.lookupOrDefault<wordRes>("patches", wordRes());

    labelList patchIDs
    (
        pbm.patchSet(patchNames).sortedToc()
    );

    labelHashSet patchSet;

    forAll(patchIDs, patchID) {
        patchSet.insert(patchIDs[patchID]);
    }

    volScalarField wallDist
    (
        patchesOutputField,
        wallDist::New(mesh, patchSet, "wall").y()
    );

    wallDist.write();

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
