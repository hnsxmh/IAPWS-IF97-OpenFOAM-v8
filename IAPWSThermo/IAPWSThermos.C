/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "psiThermo.H"
#include "rhoThermo.H"
#include "IAPWSThermo.H"
#include "pureMixture.H"

#include "forGases.H"
#include "forLiquids.H"
#include "forPolynomials.H"
#include "forTabulated.H"
#include "makeThermo.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    #define makeBasicExternalLibraryBasedThermo(Cthermo)                          \
                                                                                \
    typedef Cthermo                                                               \
        Cthermo;                                                                  \
                                                                                \
    defineTypeNameAndDebug                                                        \
    (                                                                             \
        Cthermo,                                                                  \
        0                                                                         \
    );                                                                            \
                                                                                \
    addToRunTimeSelectionTable                                                    \
    (                                                                             \
        fluidThermo,                                                           \
        Cthermo,                                                                  \
        fvMesh                                                                    \
    );     

    makeBasicExternalLibraryBasedThermo
    (
        IAPWSThermo
    );                                                                                                                                    
    // forGases(makeThermos, rhoThermo, IAPWSThermo, pureMixture);
    // forGases(makeThermos, rhoThermo, IAPWSThermo, pureMixture);
    // forLiquids(makeThermos, rhoThermo, IAPWSThermo, pureMixture);
    // forPolynomials(makeThermos, rhoThermo, IAPWSThermo, pureMixture);
    // forTabulated(makeThermos, rhoThermo, IAPWSThermo, pureMixture);
}

// ************************************************************************* //
