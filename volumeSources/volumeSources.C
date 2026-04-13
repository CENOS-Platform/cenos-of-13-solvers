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

#include "volumeSources.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(volumeSources, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        volumeSources,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::volumeSources::readCoeffs(const dictionary& dict)
{
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    UName_ =
    dict.lookupOrDefault<word>
    (
        "U",
        IOobject::groupName("U", phaseName_)
    );
    FName_ = dict.lookupOrDefault<word>("F", word::null);
    qName_ = dict.lookupOrDefault<word>("q", word::null);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::volumeSources::volumeSources
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    FName_(word::null),
    qName_(word::null),
    rho_("rho", dimDensity, dict),
    fPtr_(nullptr),
    qPtr_(nullptr)

{
    readCoeffs(coeffs(dict));
    if (FName_ != word::null)
    {
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

    if (qName_ != word::null)
    {
        qPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    qName_,
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::volumeSources::addSupFields() const
{
    wordList fields;

    if (fPtr_.valid()) fields.append(UName_);
    if (qPtr_.valid())  {
        const basicThermo& thermo =
        mesh().lookupObject<basicThermo>(physicalProperties::typeName);
        fields.append(thermo.he().name());

    }   
    return fields;
}


void Foam::fv::volumeSources::addSup
(
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    eqn += *fPtr_/rho_;
}

void Foam::fv::volumeSources::addSup
(
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    eqn += *qPtr_;
}


void Foam::fv::volumeSources::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    eqn += *fPtr_;
}


void Foam::fv::volumeSources::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    eqn += alpha * (*fPtr_);
}


bool Foam::fv::volumeSources::movePoints()
{
    return true;
}


void Foam::fv::volumeSources::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::volumeSources::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::volumeSources::distribute(const polyDistributionMap&)
{}


bool Foam::fv::volumeSources::read(const dictionary& dict)
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
