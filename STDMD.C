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

Foam::label Foam::functionObjects::STDMD::nComponents(const word &fieldName) const
{
    label nComps = 0;
    bool processed = false;

    processed = processed || nComponents<scalar>(fieldName, nComps);
    processed = processed || nComponents<vector>(fieldName, nComps);
    processed = processed || nComponents<tensor>(fieldName, nComps);

    if (!processed)
    {
        FatalErrorInFunction
            << "Unknown type of input field during initialisation: "
            << fieldName << nl
            << exit(FatalError);
    }

    return nComps;
}

// Create snapshots
void Foam::functionObjects::STDMD::snapshot()
{
    bool processed = false;
    // Check the type of the field
    processed = processed || getSnapshot<scalar>();
    processed = processed || getSnapshot<vector>();
    processed = processed || getSnapshot<tensor>();

    if (!processed)
    {
        FatalErrorInFunction
            << "    functionObjects::" << type() << " " << name() << ":"
            << " cannot find required input field during *[snapshot loading]*: "
            << fieldName_ << nl
            << "    Do you execute required functionObjects"
            << " before executing SDMD, e.g. mapFields?"
            << exit(FatalError);
    }
}

void Foam::functionObjects::STDMD::initialize()
{
    bool processed = false;

    // Check the type of the field
    processed = processed || initializeSnap<scalar>();
    processed = processed || initializeSnap<vector>();
    processed = processed || initializeSnap<tensor>();

    if (!processed)
    {
        FatalErrorInFunction
            << "    functionObjects::" << type() << " " << name() << ":"
            << " cannot find required input field during *[snapshot loading]*: "
            << fieldName_ << nl
            << exit(FatalError);
    }
}

// Classical Gram-Schmidt reorthonormalization
Foam::RectangularMatrix<double_t>
Foam::functionObjects::STDMD::GSOrthonormalize(RMatrix &x, string filePath) const
{
    RMatrix ex_ = x;
    for (int i = 0; i < nGram_; i++)
    {
        RMatrix dx_ = vectorMatrixMulti(ex_, filePath, QxCol_);
        RMatrix cx_ = vectorMatrixMulti(dx_, filePath, QyCol_);
        ex_ = ex_ - cx_;
    }
    return ex_;
}

// Expand Q, G and A based on x
void Foam::functionObjects::STDMD::expandQx(const RMatrix &ex_, const scalar exNorm_)
{
    // Qx.setSize(Qx.m(), Qx.n() + 1);
    Gx.setSize(Gx.m() + 1, Gx.n() + 1);
    A.setSize(A.m(), A.n() + 1);

    // MatrixBlock<RMatrix> qxExpand(Qx, Qx.m(), 1, 0, Qx.n() - 1);
    RMatrix qxExpand = ex_ / exNorm_;
    addCol(qxExpand, filePathQx);
    // qxExpand = ex_ / exNorm_;
}

// Expand Q, G and A based on y
void Foam::functionObjects::STDMD::expandQy(const RMatrix &ey_, const scalar eyNorm_)
{
    // Qy.setSize(Qy.m(), Qy.n() + 1);
    Gy.setSize(Gy.m() + 1, Gy.n() + 1);
    A.setSize(A.m() + 1, A.n());

    RMatrix qyExpand = ey_ / eyNorm_;
    addCol(qyExpand, filePathQy);
    // MatrixBlock<RMatrix> qyExpand(Qy, Qy.m(), 1, 0, Qy.n() - 1);
    // qyExpand = ey_ / eyNorm_;
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

// Perform Q^T*x or Q*x according to the cols of vector x
Foam::RectangularMatrix<double_t>
Foam::functionObjects::STDMD::vectorMatrixMulti(const RMatrix &x, string filePath, const label QCol) const
{
    // Assume that the size of matrix Qx is: m*n
    label xSize = x.m();
    RMatrix dx_;

    if (xSize == nSnap_)
    {
        dx_ = RMatrix(QCol, 1, Zero);
    }
    else
    {
        dx_ = RMatrix(nSnap_, 1, Zero);
    }

    std::fstream fileWrite;
    fileWrite.open(filePath, std::ios::in);

    // Perform the Q^t*x;
    string s;
    int i = 0; // Colunm counter
    while (getline(fileWrite, s))
    {
        std::istringstream sa(s);
        int j = 0; // Row counter
        while (!sa.eof())
        {
            double_t Q_ij;
            sa >> Q_ij;
            if (xSize == nSnap_) // Q*x i.e. (m*n*n*1)
            {
                dx_(i, 0) += Q_ij * x(j, 0);
            }
            else // Q^T*x i.e. (n*m*m*1)
            {
                dx_(j, 0) += Q_ij * x(i, 0);
            }
            j++;
        }
        if (j != nSnap_)
        {
            FatalErrorInFunction
                << " functionObjects::" << type() << " " << name() << ":"
                << " The rows of stored file Q_ is not consist with the number of element in a snapshot:"
                << " The rows read by fstream: " << j << nl
                << " The number of element in a snapshot: " << nSnap_ << nl
                << exit(FatalError);
        }
        s.clear();
        i++;
    }
    // Info << "Col of file Q is: " << i << endl;
    // Info << "The record of Q is: " << QxCol_ <<endl;
    fileWrite.close();
    return dx_;
}

//  Add vector to Qx and Qy
void Foam::functionObjects::STDMD::addCol(const RMatrix &x, string filePath) const
{
    ofstream fileWrite;
    fileWrite.open(filePath, std::ios::app);
    label endOF = x.m() - 1;
    forAll(x, i)
    {
        fileWrite << x(i, 0);
        if (i != endOF)
        {
            fileWrite << " ";
        }
    }
    fileWrite << "\n";
    fileWrite.close();
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
      QxCol_(1),
      QyCol_(1),
      A(),
      Gx(),
      Gy(),
      x_(),
      y_(),
      fieldName_(dict.lookupOrDefault<word>("field", "U")),
      nGram_(5),
      outputDir_()
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
      QxCol_(1),
      QyCol_(1),
      A(),
      Gx(),
      Gy(),
      x_(),
      y_(),
      fieldName_(dict.lookupOrDefault<word>("field", "U")),
      nGram_(5),
      outputDir_()
{
    // Read settings from dictionary files
    read(dict);
}

// * * * * * * * * * Destructor  * * * * * * * * * * * * * //
Foam::functionObjects::STDMD::~STDMD()
{
}

// * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::STDMD::read(const dictionary &dict)
{
    fvMeshFunctionObject::read(dict);
    Log << "The type of post-processing:" << type() << endl;
    Info << "The selected field: " << fieldName_ << endl;

    // Read locations of two points to limit the aera
    dict.lookup("pointLocations") >> pointLocation_;

    return true;
}

// Write A, Gx, Gy, Qx, Qy to files
bool Foam::functionObjects::STDMD::write()
{

    if (step_ > 2)
    {
        Info << "Writing Matrix to postProcessing" << endl;

        // fileName outputDir = mesh_.time().path() / ".." / "postProcessing" / "SDMD";
        writeMatrix(outputDir_, A, "A");
        // writeMatrix(outputDir_, Qx, "Qx");
        // writeMatrix(outputDir_, Qy, "Qy");
        writeMatrix(outputDir_, Gx, "Gx");
        writeMatrix(outputDir_, Gy, "Gy");
    }
    return true;
}

bool Foam::functionObjects::STDMD::execute()
{
    //  Output relevant information about SDMD
    Log << "Execute Function is running" << endl;
    Log << type() << " " << name() << " execute:" << endl;

    // SDMD is processed by the master processor
    Log << "Execution index:" << step_ << endl;

    if (step_ == 0)
    {
        initialize();
        snapshot();
        writeMatrix(outputDir_, y_, "snapshots1");
    }
    else
    {
        if (Pstream::master())
        {
            x_ = y_;
        }
        snapshot();

        // Process the First Iterate
        if (step_ == 1)
        {
            Info << "=======Step 1=======" << endl;
            if (Pstream::master())
            {
                scalar normX_ = L2norm(x_);
                scalar normY_ = L2norm(y_);

                RMatrix QxFirstCol = x_ / normX_;
                RMatrix QyFirstCol = y_ / normY_;

                addCol(QxFirstCol, filePathQx);
                addCol(QyFirstCol, filePathQy);

                Gx(0, 0) = normX_ * normX_;
                Gy(0, 0) = normY_ * normY_;
                A(0, 0) = normX_ * normY_;
            }
        }

        if (Pstream::master())
        {
            // Algorithm step 1
            // i.e. Gram-Schmidt reothonormalization
            RMatrix ex_ = GSOrthonormalize(x_, filePathQx);
            scalar normEx_ = L2norm(ex_);

            RMatrix ey_ = GSOrthonormalize(y_, filePathQy);
            scalar normEy_ = L2norm(ey_);

            // Algorithm step 2
            // Check basis for x_ and expand, if necessary
            scalar normX_ = L2norm(x_);
            scalar normY_ = L2norm(y_);

            if (normEx_ / normX_ > __DBL_EPSILON__)
            {
                Info << "Expand the Qx" << endl;
                expandQx(ex_, normEx_);
                QxCol_++;
            }

            if (normEy_ / normY_ > __DBL_EPSILON__)
            {
                Info << "Expand the Qy" << endl;
                expandQy(ey_, normEy_);
                QyCol_++;
            }
            // Algorithm step 3
            // Check if POD compression is need
            // This step is not considered for now

            // Algorithm step 4
            // Calculate xtilde and ytilde
            Info << "Caculate tilde" << endl;
            RMatrix xtilde_ = vectorMatrixMulti(x_, filePathQx, QxCol_);
            RMatrix ytilde_ = vectorMatrixMulti(y_, filePathQy, QyCol_);
            Info << "The size of xtilde is: " << xtilde_.m()<<endl;
            Info << "The size of ytilde is: " << ytilde_.m()<<endl;
            // RMatrix xtilde_ = calcTilde(Qx, x_);
            // RMatrix ytilde_ = calcTilde(Qy, y_);

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
