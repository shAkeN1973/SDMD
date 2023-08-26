/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"

// Libraries declared here are transient and may affect compliation
#include "fvMeshFunctionObject.H"
#include "STDMD.H"

// Initialize the first snapshot and matrix utilized in SDMD
template <class Type>
bool Foam::functionObjects::STDMD::initializeSnap()
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    
    if (mesh_.foundObject<volFieldType>(fieldName_))
    {
        // Get velocity field reference
        const volFieldType &field = lookupObject<volFieldType>(fieldName_);

        const label nComps = nComponents(fieldName_);
        Log << "Components is: " << nComps << endl;

        nSnap_ = nComps * mesh_.nCells();
        reduce(nSnap_, sumOp<label>());

        nCells_ = mesh_.nCells();
        reduce(nCells_, sumOp<label>());

        outputDir_ = mesh_.time().path() / ".." / "postProcessing" / "SDMD";

        if (Pstream::parRun())
        {
            // Gather the central point together
            List<pointField> listMeshCentre(Pstream::nProcs());
            listMeshCentre[Pstream::myProcNo()] = field.mesh().C();
            Pstream::gatherList(listMeshCentre);
            centralPoint_ = pointField(nCells_, Zero);

            if (Pstream::master())
            {
                label iMeshBlockLen = 0; // Length of each mesh block
                // Initialize the coordinates of mesh points in centralPoint
                for (int iProc = 0; iProc < Pstream::nProcs(); iProc++)
                {
                    pointField &ct = listMeshCentre[iProc];
                    for (int iCell = 0; iCell < listMeshCentre[iProc].size(); iCell++)
                    {
                        centralPoint_[iCell + iMeshBlockLen] = ct[iCell];
                    }
                    iMeshBlockLen += listMeshCentre[iProc].size();
                }

                // Check if process aera is limited
                if (pointLocation_.size())
                {
                    nCells_ = 0; // Resize the number of all meshs
                    double xMax = pointLocation_[0].x();
                    double yMax = pointLocation_[0].y();
                    double zMax = pointLocation_[0].z();

                    double xMin = pointLocation_[1].x();
                    double yMin = pointLocation_[1].y();
                    double zMin = pointLocation_[1].z();

                    pointField tempCentralPoint;
                    forAll(centralPoint_, iPoint)
                    {
                        point &ct = centralPoint_[iPoint];
                        if ((ct.x() <= xMax && ct.x() >= xMin) && (ct.y() <= yMax && ct.y() >= yMin) && (ct.z() <= zMax && ct.z() >= zMin))
                        {
                            nCells_++;
                            pointIndexList_.append(iPoint);
                            tempCentralPoint.append(ct);
                        }
                    }
                    centralPoint_ = tempCentralPoint;
                    nSnap_ = nComps * nCells_;
                }

                Log << "Cells of total mesh cells: " << nCells_ << endl;
                Log << "Number of elements in a snapshot: " << nSnap_ << endl;

                if (nSnap_ > 0)
                {
                    x_ = RMatrix(nSnap_, 1, Zero);
                    y_ = RMatrix(nSnap_, 1, Zero);
                    Qx = RMatrix(nSnap_, 1, Zero);
                    Qy = RMatrix(nSnap_, 1, Zero);
                    Gx = RMatrix(1, 1, Zero);
                    Gy = RMatrix(1, 1, Zero);
                    A = RMatrix(1, 1, Zero);
                }
                else
                {
                    x_ = RMatrix(1, 1, Zero);
                    y_ = RMatrix(1, 1, Zero);
                    Qx = RMatrix(1, 1, Zero);
                    Qy = RMatrix(1, 1, Zero);
                    Gx = RMatrix(1, 1, Zero);
                    Gy = RMatrix(1, 1, Zero);
                    A = RMatrix(1, 1, Zero);
                }

                // Write the coordinates of central points
                mkDir(outputDir_);
                OFstream osCoordinate(
                    outputDir_ / "coordinate.raw",
                    IOstream::ASCII,
                    IOstream::currentVersion,
                    IOstream::UNCOMPRESSED);

                forAll(centralPoint_, elementi)
                {
                    point &pt = centralPoint_[elementi];
                    osCoordinate << pt.x() << " "
                                 << pt.y() << " "
                                 << pt.z() << endl;
                }
            }
        }
        return true;
    }
    return false;
}

// Get snapshots from data field
template <class Type>
bool Foam::functionObjects::STDMD::getSnapshot()
{
    // typedef GeometricField<vector, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef Field<Type> fieldType;

    if (mesh_.foundObject<volFieldType>(fieldName_))
    {
        const label nComps = nComponents(fieldName_);

        // Get size of Vector field
        const volFieldType &field = lookupObject<volFieldType>(fieldName_);
        // const label nField = field.size();

        // Get the number of element  of mesh and one snapshot in original aera
        label nSnapTotal = nComps * mesh_.nCells();
        reduce(nSnapTotal, sumOp<label>());

        label nCellsTotal = mesh_.nCells();
        reduce(nCellsTotal, sumOp<label>());

        // Parallel process:
        if (Pstream::parRun())
        {
            // Gather list from openfoam field
            List<fieldType> listVectorField(Pstream::nProcs());
            listVectorField[Pstream::myProcNo()] = field.internalField();
            Pstream::gatherList(listVectorField);

            Info << "Elements in original aera: " << nSnapTotal << endl;
            if (Pstream::master())
            {
                //  Gather vectorField from vectorFieldList
                RMatrix tempY(nSnapTotal, 1, Zero);

                // Assign velocity to snapshots matrix x_
                for (direction dir = 0; dir < nComps; dir++)
                {
                    label mStart_ = 0;
                    for (int iProc = 0; iProc < Pstream::nProcs(); iProc++)
                    {
                        fieldType &UField = listVectorField[iProc];
                        const label iProcFieldSize_ = UField.size();
                        MatrixBlock<RMatrix> v(tempY, iProcFieldSize_, 1, mStart_ + dir * nCellsTotal, 0);
                        v = UField.component(dir);
                        mStart_ += UField.size();
                    }
                }
                // Limit aera process:
                if (pointLocation_.size())
                {
                    for (direction dir = 0; dir < nComps; dir++)
                    {
                        for (int i = 0; i < pointIndexList_.size(); i++)
                        {
                            label index = pointIndexList_[i];
                            y_(i + dir * nCells_, 0) = tempY(index + dir * nCellsTotal, 0);
                        }
                    }
                }
                else 
                {
                    y_ = tempY;
                }
                Info << "size point index list: " << pointIndexList_.size() << endl;
                Info << "size of resized mesh: " << nCells_<<endl;
            }
            return true;
        }
    }

    return false;
}

// Returns the number of components in the seleted field
template <class Type>
bool Foam::functionObjects::STDMD::nComponents(const word &fieldName, label &nComps) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    if (mesh_.foundObject<volFieldType>(fieldName))
    {
        nComps = pTraits<typename volFieldType::value_type>::nComponents;
        return true;
    }

    return false;
}

// ************************************************************************* //