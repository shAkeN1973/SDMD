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

Description
    A Streaming algorithm used for the function objects

    functions: 
        - Performing the DMD to the new snapshots

Class
    Foam::functionObjects::SDMD.H
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
:
    public fvMeshFunctionObject,
    public logFiles
{
    typedef RectangularMatrix<double_t>  RMatrix;
    typedef SquareMatrix<double_t> SMatrix;

// Private Data

    // Protected Data

    // Number of elements in a snapshot
    label nSnap_;
    
    // Number of the snapshots
    label snapNum_;

    // Maximum allowable rank of DMD operator
    label maxRank_;

    // Number of pairs of snapshots processed
    label step_;

    // Q_x, whose columns form an orthonormal basis for the image of X
    // Eq: (1)
    RMatrix Qx;

    // Q_y, whose columns form an orthonormal basis for the image of Y
    // Eq: (1)
    RMatrix Qy;

    // A, New matrix defined
    SMatrix A;

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

// Private Member Functions

    // Number of the field components, e.g. pressure filed 1, Velocity field 3
    label nComponents(const word& fieldName);

    // Return L2-norm of a given column vector
    scalar L2norm(const RMatrix& z) const;

    // Get snapshot form the Volum field
    bool getSnapshot();

    // Initialise the first snapshots
    void initialize();

    // Classical Gram-Schmidt reorthonormalization 
    RMatrix GSOrthonormalize(RMatrix x,RMatrix Q) const;

public:
    // Runtime type information
    TypeName("STDMD");

    // Constructors
    STDMD
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    STDMD
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

    // Destructor
    virtual ~STDMD();

    // Member Functions

    // Read the dict data
    virtual bool read(const dictionary&);
    virtual bool execute();
    virtual bool write();
    virtual void writeFileHeader(const label i);

};

}   // End namespace functionObjects
}   // End namespace Foam

#endif



