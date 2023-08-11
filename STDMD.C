/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    SDMD

Description
    Streaming DMD for openFOAM v7
\*---------------------------------------------------------------------------*/

#include "STDMD.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //
namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(STDMD, 0);

        addToRunTimeSelectionTable(functionObject, STDMD, dictionary);
    }
}
// * * * * * * * * * Protected Member Functions* * * * * * * * * * * * * //

// * * * * * * * * * Private Member Functions* * * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::STDMD::L2norm(const RMatrix &z) const
{
#ifdef FULLDEBUG
    // Check if the given RectangularMatrix is effectively a column vector
    if (z.n() != 1)
    {
        FatalErrorInFunction
            << "Input matrix is not a column vector."
            << exit(FatalError);
    }
#endif

    label mRows_ = z.m();
    // EInfo << "Rows of input matrix: " << mRows_ << endl;

    scalar result = Zero;

    // L2 norm i.e. : sqrt(sum{z(i)^2})
    for (label i = 0; i < mRows_; ++i)
    {
        result += magSqr(z(i, 0));
    }

    // reduce(result, sumOp<scalar>());
    //  Heuristic addition to avoid very small or zero norm
    return max(SMALL, Foam::sqrt(result));
}

void Foam::functionObjects::STDMD::initialize()
{
    const label nComps = nComponents(fieldName_);
    Log << "Components is: " << nComps << endl;

    nSnap_ = nComps * mesh_.nCells();
    reduce(nSnap_, sumOp<label>());

    nCells_ = mesh_.nCells();
    reduce(nCells_, sumOp<label>());

    Log << "Number of elements in one Snapshots:" << nSnap_ << endl;

    // Get velocity field reference
    const volVectorField &field = lookupObject<volVectorField>(fieldName_);

    if (Pstream::parRun())
    {
        // Gather the central point together
        List<pointField> listMeshCentre(Pstream::nProcs());
        listMeshCentre[Pstream::myProcNo()] = field.mesh().C();
        Pstream::gatherList(listMeshCentre);

        Log << "Cells of total mesh cells: " << nCells_ << endl;

        centralPoint_ = pointField(nCells_, Zero);

        if (Pstream::master())
        {
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

            // Write the coordinates of central points
            fileName outputDir = mesh_.time().path() / ".." / "postProcessing" / "SDMD";

            mkDir(outputDir);
            OFstream osCoordinate(
                outputDir / "coordinate.raw",
                IOstream::ASCII,
                IOstream::currentVersion,
                IOstream::UNCOMPRESSED);

            // Header
            osCoordinate << "Grid center point coordinates" << nl
                         << "x "
                         << "y "
                         << "z" << endl;

            forAll(centralPoint_, elementi)
            {
                point &pt = centralPoint_[elementi];
                osCoordinate << pt.x() << " "
                             << pt.y() << " "
                             << pt.z() << endl;
            }
        }
    }
}

// Returns the number of components in the seleted field
Foam::label Foam::functionObjects::STDMD::nComponents(const word &fieldName)
{
    label nComps;
    typedef GeometricField<vector, fvPatchField, volMesh> VolFieldType;
    if (mesh_.foundObject<VolFieldType>(fieldName))
    {
        nComps = pTraits<typename VolFieldType::value_type>::nComponents;
        return nComps;
    }
    return 0;
}

// Get snapshots from data field
bool Foam::functionObjects::STDMD::getSnapshot()
{

    const label nComps = nComponents(fieldName_);
    // typedef GeometricField<vector, fvPatchField, volMesh> VolFieldType;

    // Get size of Vector field
    const volVectorField &field = lookupObject<volVectorField>(fieldName_);
    const label nField = field.size();

    // Parallel process:
    if (Pstream::parRun())
    {
        // Gather list from openfoam field
        List<vectorField> listVectorField(Pstream::nProcs());
        listVectorField[Pstream::myProcNo()] = field.internalField();
        Pstream::gatherList(listVectorField);

        if (Pstream::master())
        {
            //  Gather vectorField from vectorFieldList
            const label nComps_ = nComponents(fieldName_);

            // Assign velocity to snapshots matrix x_
            for (direction dir = 0; dir < nComps_; dir++)
            {
                label mStart_ = 0;
                for (int iProc = 0; iProc < Pstream::nProcs(); iProc++)
                {
                    vectorField &UField = listVectorField[iProc];
                    const label iProcFieldSize_ = UField.size();
                    MatrixBlock<RMatrix> v(y_, iProcFieldSize_, 1, mStart_ + dir * nCells_, 0);
                    v = UField.component(dir);
                    mStart_ += UField.size();
                }
            }
        }
        return true;
    }
    else // Non-Parallel processing
    {
        // Assignment the U vector to matrix x_
        for (direction dir = 0; dir < 3; dir++)
        {
            MatrixBlock<RMatrix> v(y_, nField, 1, 0 + dir * nField, 0);
            v = field.component(dir);
        }

        fileName outputDir = mesh_.time().path() / "Postprocessing";
        mkDir(outputDir);

        OFstream osNonParallel(
            outputDir / "NonParallell_U_Mesh.raw",
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED);

        pointField meshCentr = field.mesh().C();
        vectorField Umesh = field.internalField();

        // Write the central coordinate of each mesh cell as well as velocity component
        for (int i = 0; i < field.mesh().nCells(); i++)
        {
            point &ct = meshCentr[i];
            vector &Uv = Umesh[i];
            osNonParallel << ct.x() << " " << ct.y() << " " << ct.z() << " "
                          << Uv.x() << " " << Uv.y() << " " << Uv.z() << " " << nl;
        }
        return true;
    }
    return false;
}

// Classical Gram-Schmidt reorthonormalization
Foam::RectangularMatrix<double_t>
Foam::functionObjects::STDMD::GSOrthonormalize(RMatrix &x, RMatrix &Q) const
{

    RMatrix ex_ = x;

    for (int i = 0; i < nGram_; i++)
    {
        // dx = Q^T * ex;
        RMatrix dx_(Q.n(), 1, Zero);
        for (int col_ = 0; col_ < Q.n(); col_++)
        {
            for (int row_ = 0; row_ < Q.m(); row_++)
            {
                dx_(col_, 0) += Q(row_, col_) * ex_(row_, 0);
            }
        }
        ex_ = ex_ - Q * dx_;
    }
    return ex_;
}

// Expand Q, G and A based on x
void Foam::functionObjects::STDMD::expandQx(const RMatrix &ex_, const scalar exNorm_)
{
    Qx.setSize(Qx.m(), Qx.n() + 1);
    Gx.setSize(Gx.m() + 1, Gx.n() + 1);
    A.setSize(A.m(), A.n() + 1);

    MatrixBlock<RMatrix> qxExpand(Qx, Qx.m(), 1, 0, Qx.n() - 1);
    qxExpand = ex_ / exNorm_;
}

// Expand Q, G and A based on y
void Foam::functionObjects::STDMD::expandQy(const RMatrix &ey_, const scalar eyNorm_)
{
    Qy.setSize(Qy.m(), Qy.n() + 1);
    Gy.setSize(Gy.m() + 1, Gy.n() + 1);
    A.setSize(A.m() + 1, A.n());

    MatrixBlock<RMatrix> qyExpand(Qy, Qy.m(), 1, 0, Qy.n() - 1);
    qyExpand = ey_ / eyNorm_;
}

// Calculate xtilde
Foam::RectangularMatrix<double_t>
Foam::functionObjects::STDMD::calcTilde(RMatrix &Q, RMatrix &x) const
{
    RMatrix xtilde_(Q.n(), 1, Zero);
    for (int col_ = 0; col_ < Q.n(); col_++)
    {
        for (int row_ = 0; row_ < Q.m(); row_++)
        {
            xtilde_(col_, 0) += Q(row_, col_) * x(row_, 0);
        }
    }
    return xtilde_;
}

// Return the transpose of the matrix
Foam::RectangularMatrix<double_t>
Foam::functionObjects::STDMD::transpose(const RMatrix &A) const
{
    RMatrix T(A.n(), A.m(), Zero);
    for (int row_ = 0; row_ < A.m(); row_++)
    {
        for (int col_ = 0; col_ < A.n(); col_++)
        {
            T(col_, row_) = A(row_, col_);
        }
    }
    return T;
}

// Write the specified matrix to postProcessing/SDMD
void Foam::functionObjects::STDMD::writeMatrix(
    const fileName &outputDir,
    const RMatrix &A,
    const fileName &matrixName) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }
    OFstream os(
        outputDir / matrixName + ".raw",
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED);

    for (int row_ = 0; row_ < A.m(); row_++)
    {
        for (int col_ = 0; col_ < A.n(); col_++)
        {
            os << A(row_, col_) << " ";
        }
        os << endl;
    }
}

// * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::functionObjects::STDMD::STDMD(
    const word &name,
    const Time &runTime,
    const dictionary &dict)
    : fvMeshFunctionObject(name, runTime, dict),
      logFiles(obr_, name),
      nSnap_(),
      snapNum_(),
      maxRank_(),
      step_(0),
      Qx(),
      Qy(),
      A(),
      Gx(),
      Gy(),
      x_(),
      y_(),
      fieldName_(dict.lookupOrDefault<word>("field", "U")),
      nGram_(5)
{
    // Read settings from dictionary files
    read(dict);
}

Foam::functionObjects::STDMD::STDMD(
    const word &name,
    const objectRegistry &obr,
    const dictionary &dict)
    : fvMeshFunctionObject(name, obr, dict),
      logFiles(obr_, name),
      nSnap_(),
      snapNum_(),
      maxRank_(),
      step_(0),
      Qx(),
      Qy(),
      A(),
      Gx(),
      Gy(),
      x_(),
      y_(),
      fieldName_(dict.lookupOrDefault<word>("field", "U")),
      nGram_(5)
{
    // Read settings from dictionary files
    read(dict);
}
// * * * * * * * * * Destructor  * * * * * * * * * * * * * //

Foam::functionObjects::STDMD::~STDMD()
{
}

// * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// mesh_ : refrence to fvMesh

bool Foam::functionObjects::STDMD::read(const dictionary &dict)
{
    fvMeshFunctionObject::read(dict);

    Log << "This is type:" << type() << nl
        << "this is name:" << name() << ":" << nl;

    return true;
}

// Write A, Gx, Gy, Qx, Qy to files
bool Foam::functionObjects::STDMD::write()
{

    if (step_ > 2)
    {
        Info << "Writing Matrix to postProcessing" << endl;

        fileName outputDir = mesh_.time().path() / ".." / "postProcessing" / "SDMD";

        writeMatrix(outputDir, A, "A");
        writeMatrix(outputDir, Qx, "Qx");
        writeMatrix(outputDir, Qy, "Qy");
        writeMatrix(outputDir, Gx, "Gx");
        writeMatrix(outputDir, Gy, "Gy");
    }

    return true;
}

bool Foam::functionObjects::STDMD::execute()
{
    // Log information
    Log << "execute Function is running" << endl;
    Log << type() << " " << name() << " execute:" << endl;

    // SDMD is processed by the master processor

    Log << "Execution index:" << step_ << endl;

    if (step_ == 0)
    {
        initialize();
        getSnapshot();
    }
    else
    {
        if (Pstream::master())
        {
            x_ = y_;
        }
        getSnapshot();

        // Process the First Iterate
        if (step_ == 1)
        {
            Info << "=======Step 1=======" << endl;
            if (Pstream::master())
            {
                scalar normX_ = L2norm(x_);
                scalar normY_ = L2norm(y_);

                Qx = x_ / normX_;
                Qy = y_ / normY_;
                Gx(0, 0) = normX_ * normX_;
                Gy(0, 0) = normY_ * normY_;
                A(0, 0) = normX_ * normY_;
            }
        }

        if (Pstream::master())
        {
            // Algorithm step 1
            // i.e. Gram-Schmidt reothonormalization
            RMatrix ex_ = GSOrthonormalize(x_, Qx);
            scalar normEx_ = L2norm(ex_);

            RMatrix ey_ = GSOrthonormalize(y_, Qy);
            scalar normEy_ = L2norm(ey_);

            // Algorithm step 2
            // Check basis for x_ and expand, if necessary
            scalar normX_ = L2norm(x_);
            scalar normY_ = L2norm(y_);

            if (normEx_ / normX_ > __DBL_EPSILON__)
            {
                expandQx(ex_, normEx_);
            }

            if (normEy_ / normY_ > __DBL_EPSILON__)
            {
                expandQy(ey_, normEy_);
            }
            // Algorithm step 3
            // Check if POD compression is need
            // This step is not considered for now

            // Algorithm step 4
            // Calculate xtilde and ytilde
            Info << "Caculate tilde" << endl;
            RMatrix xtilde_ = calcTilde(Qx, x_);
            RMatrix ytilde_ = calcTilde(Qy, y_);

            // Update A and Gx,Gy
            Info << "Update A and Gx, Gy" << endl;
            A = A + ytilde_ * transpose(xtilde_);
            Gx = Gx + xtilde_ * transpose(xtilde_);
            Gy = Gy + ytilde_ * transpose(ytilde_);
        }
    }
    step_++;
    return true;
}

void Foam::functionObjects::STDMD::writeFileHeader(const label i)
{
    Log << "write File Header is running";
}

// ************************************************************************* //
