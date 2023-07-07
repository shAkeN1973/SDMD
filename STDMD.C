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
    FOtest

Description
    A test file used for the function objects
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
    defineTypeNameAndDebug(STDMD,0);

    addToRunTimeSelectionTable(functionObject,STDMD,dictionary);
}
}
// * * * * * * * * * Protected Member Functions* * * * * * * * * * * * * //



// * * * * * * * * * Private Member Functions* * * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::STDMD::L2norm(const RMatrix& z) const
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
    scalar result = Zero;

    // L2 norm i.e. : sqrt(sum{z(i)^2})
    for(label i = 0; i < mRows_; ++i){
        result += magSqr(z(i,0));
    }
    reduce(result, sumOp<scalar>());
    // Heuristic addition to avoid very small or zero norm
    return max(SMALL, Foam::sqrt(result));
}


void Foam::functionObjects::STDMD::initialize(){
    Log<<"Initialize is runnging"<<endl;
    const label nComps = nComponents(fieldName_);

    Log<<"Components is: "<<nComps<<endl;

    nSnap_ = nComps*mesh_.nCells();

    if(nSnap_>0){
        z_ = RMatrix(nSnap_,1,Zero);
    }
    else{
        z_ = RMatrix(1,1,Zero);
    }
    Log<<"Number of elements in one Snapshots:"<<nSnap_<<endl;
}

// Returns the number of components in the seleted field
Foam::label Foam::functionObjects::STDMD::nComponents(const word &fieldName){
    label nComps;
    typedef GeometricField<vector, fvPatchField, volMesh> VolFieldType;
    if(mesh_.foundObject<VolFieldType>(fieldName)){
        nComps = pTraits<typename VolFieldType::value_type>::nComponents;
        return nComps;
    }
    return 0;
}

// Get snapshots from data field
bool Foam::functionObjects::STDMD::getSnapshot(){
    const label nComps = nComponents(fieldName_);
    typedef GeometricField<vector, fvPatchField, volMesh> VolFieldType;
    const VolFieldType& field = lookupObject<VolFieldType>(fieldName_);
    const label nField = field.size();
    Log<<"field Size of the U vector:"<<nField<<endl;
    direction dir;
    dir = 0;
   MatrixBlock<RMatrix> v(z_,nField,1,0,0);
   v = field.component(dir);
   Info<<"Number of rows in block v:"<<v.m()<<nl
   <<"Number of columns in block v:"<<v.n()<<endl;

   // Write files

   // Create the ourput file dictionary
   fileName outputPath_ = mesh_.time().path()/"postProcessing";
   mkDir(outputPath_);
    
   //autoPtr<OFstream> outputFilePtr;
   OFstream os
   (
        outputPath_/"123.raw",
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED 
   );

    os  << "This is header:"<<nl
    <<"content of v vector"<<nl;

    for(label i=0; i<nField; i++){
        os<< v(0,i)<<nl;
    };



    return false;
}


// * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::functionObjects::STDMD::STDMD
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name,runTime,dict),
    logFiles(obr_,name),
    nSnap_(),
    snapNum_(),
    maxRank_(),
    step_(0),
    Qx(),
    Qy(),
    A(),
    Gx(),
    Gy(),
    z_(),
    fieldName_(dict.lookupOrDefault<word>("field","U"))
{
    // Read settings from dictionary files
    read(dict);
}

Foam::functionObjects::STDMD::STDMD
(
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
)
:
fvMeshFunctionObject(name,obr,dict),
logFiles(obr_,name),
nSnap_(),
snapNum_(),
maxRank_(),
step_(0),
Qx(),
Qy(),
A(),
Gx(),
Gy(),
z_(),
fieldName_(dict.lookupOrDefault<word>("field","U"))
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


bool Foam::functionObjects::STDMD::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log<<"This is type:"<<type()<<nl
    <<"this is name:"<<name()<<":"<<nl;

    return true;
}

bool Foam::functionObjects::STDMD::write()
{
    return true;
}

bool Foam::functionObjects::STDMD::execute()
{
    // Log information 
    Log<<"execute Function is running"<< endl;
    Log << type() << " " << name() << " execute:" << endl;

    if(step_ == 0){
        initialize();
    }
    
    if(step_ > 0){
        Log<<"Execution index:"<<step_<< endl;
        getSnapshot();
    }

    step_++;
 
    return true;
}

void Foam::functionObjects::STDMD::writeFileHeader(const label i)
{
    Log<<"write File Header is running";
}


// ************************************************************************* //
