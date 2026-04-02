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

Class
    Foam::fv::volumeForce

Description
    Applies the volumetric force to the momentum equation.

Usage
    Example usage:
    \verbatim
    volumeForce
    {
        F           F;  // Name of the force field (needs to be present in case)
    }
    \endverbatim

SourceFiles
    volumeForce.C

\*---------------------------------------------------------------------------*/

#ifndef volumeForce_H
#define volumeForce_H

#include "fvModel.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class volumeForce Declaration
\*---------------------------------------------------------------------------*/

class volumeForce
:
    public fvModel
{
    // Private Data

        word phaseName_;

        //- Name of the velocity field
        word UName_;
        //- Name of the force field
        word FName_;

        dimensionedScalar rho_;
		
		autoPtr<volVectorField> fPtr_;

    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("volumeForce");


    // Constructors

        //- Construct from explicit source name and mesh
        volumeForce
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        volumeForce(const volumeForce&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Add explicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const volumeForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
