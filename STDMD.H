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
    STDMD

Class
    Foam::functionObjects::SDMD.H

Description
    A Streaming algorithm used for the function objects

    functions:
        - Performing the DMD to the new snapshots
        - Limit the aera to perform SDMD

Usage
    Minimal example by using system/controlDict.functions:

    SDMD
    {
        type            STDMD;
        functionObjectLibs ("libSTDMD.so");
        writeControl timeStep;
        writeInterval 10;
        // timeStart 1.9835;
        // timeEnd 2.9850;
        pointLoactions
        (
            (0 0 0)
            (1 1 1)
        )
        nGramSchmidt 5;
    }
\*---------------------------------------------------------------------------*/

#ifndef STDMD_H
#define STDMD_H

// include the functionObject library
#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "RectangularMatrix.H"
#include "QRMatrix.H"

namespace Foam
{
    namespace functionObjects
    {

        /*---------------------------------------------------------------------------*\
                                 Class forceCoeffs Declaration
        \*---------------------------------------------------------------------------*/

        class STDMD
            : public fvMeshFunctionObject,
              public logFiles
        {
            typedef RectangularMatrix<double_t> RMatrix;
            typedef SquareMatrix<double_t> SMatrix;

        // Private Data

            //===  MASTER CORE ===
            // Number of elements in a snapshot
            label nSnap_;

            // Number of the snapshots, i.e. the number of steps
            label snapNum_;

            // Maximum allowable rank of DMD operator
            label maxRank_;

            // Number of pairs of snapshots processed
            label step_;

            // Q_x, whose columns form an orthonormal basis for the image of X
            // Eq: (1)
          // RMatrix Qx;

            // Q_y, whose columns form an orthonormal basis for the image of Y
            // Eq: (1)
           // RMatrix Qy;

            // The number of columns of Qx
            label QxCol_;

            // The number of columns of Qy
            label QyCol_;

            // A, New matrix defined
            RMatrix A;

            // G_x
            RMatrix Gx;

            // G_y
            RMatrix Gy;

            // Snapshot pairs (x,y)
            RMatrix x_;

            // Snapshot pairs (x,y)
            RMatrix y_;

            // Name of operand field
            const word fieldName_;

            // Coordinates of central point for each cell
            pointField centralPoint_;

            // Cell number of all meshs
            label nCells_;

            // Number of times to reapply Gram-Schmidt
            const label nGram_;

            //  Output path
            fileName outputDir_;

            // Locations of points limiting the aera to perform SDMD
            pointField pointLocation_;

            // Path of Qx
            string filePathQx;

            // Path of Qy
            string filePathQy;


            // === ON CORE ===
            // List of index point in limit aera
            List<label> pointIndexListSlave_;
            
         // Private Member Functions

            // Number of the field components, e.g. pressure filed 1, Velocity field 3
            // Return number of components of the base type of a given field
            label nComponents(const word &fieldName) const;

            // Get the number of components of the base type of a given field
            template <class Type>
            bool nComponents(const word &fieldName, label &nComps) const;

            // Return L2-norm of a given column vector
            scalar L2norm(const RMatrix &z) const;

            // Create snapshots
            void snapshot();

            // Get snapshot form the Volum field
            template <class Type>
            bool getSnapshot();

            // Initialise the first snapshots
            void initialize();

            // Initialise the first snapshots
            template <class Type>
            bool initializeSnap();

            // Classical Gram-Schmidt reorthonormalization
            RMatrix GSOrthonormalize(RMatrix &x, string filePath) const;

            // Expand Qx, A, Gx
            void expandQx(const RMatrix &ex_, const scalar exNorm_);

            // Expand Qy, A, Gy
            void expandQy(const RMatrix &ey_, const scalar eyNorm_);

            // Calculate xtilde
            RMatrix calcTilde(RMatrix &Q, RMatrix &x) const;

            // Return the transpose of the matrix
            RMatrix transpose(const RMatrix &A) const;

            // Write the specificed matrix to postProcessing/SDMD
            void writeMatrix(const fileName &outputDir, const RMatrix &A, const fileName &matrixName) const;

            // Perform Q^T*x or Q*x according to the cols of vector x
            RMatrix vectorMatrixMulti(const RMatrix &x, string filePath, const label QCol) const;

            // Add  vector to Qx and Qy
            void addCol(const RMatrix &x, string filePath) const;

        public:
            // Runtime type information
            TypeName("STDMD");

            // Constructors
            STDMD(
                const word &name,
                const Time &runTime,
                const dictionary &dict);

            STDMD(
                const word &name,
                const objectRegistry &obr,
                const dictionary &dict);

            // Destructor
            virtual ~STDMD();

            // Member Functions

            // Read the dict data
            virtual bool read(const dictionary &);
            virtual bool execute();
            virtual bool write();
            virtual void writeFileHeader(const label i);
        };

    } // End namespace functionObjects
} // End namespace Foam

#ifdef NoRepository
#include "SDMDTemplate.C"
#endif

#endif
