/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "volumeForce.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(volumeForce, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        volumeForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeForce::readCoeffs(const dictionary& dict)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    UName_ =
    dict.lookupOrDefault<word>
    (
        "U",
        IOobject::groupName("U", phaseName_)
    );
    FName_ = dict.lookupOrDefault<word>("F", word::null);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumeForce::volumeForce
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    FName_(word::null),
    rho_("rho", dimDensity, dict)
    fPtr_(nullptr),

{
    readCoeffs(coeffs(dict));
    fPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                FName_,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::volumeForce::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::volumeForce::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    eqn += *fPtr_/rho_;
}


void Foam::fv::volumeForce::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    eqn += *fPtr_;
}


void Foam::fv::volumeForce::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    eqn += alpha * (*fPtr_);
}


bool Foam::fv::volumeForce::movePoints()
{
    return true;
}


void Foam::fv::volumeForce::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::volumeForce::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::volumeForce::distribute(const polyDistributionMap&)
{}


bool Foam::fv::volumeForce::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
