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

Class
    Foam::IAPWSThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    IAPWSThermo.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSThermo_H
#define IAPWSThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IAPWSThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class IAPWSThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private data

    //- DensityField
        volScalarField rho_;

        volScalarField drhodh_;

        volScalarField psi_;

    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

public:

    //- Runtime type information
    TypeName("IAPWSThermo");


    // Constructors

        //- Construct from mesh and phase name
        IAPWSThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Disallow default bitwise copy construction
        IAPWSThermo(const IAPWSThermo<BasicPsiThermo, MixtureType>&) = delete;


    //- Destructor
    virtual ~IAPWSThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const IAPWSThermo<BasicPsiThermo, MixtureType>&
        ) = delete;

    // Fields derived from thermodynamic state variables

        //- Enthalpy for cell-set [J/kg]
        virtual tmp<scalarField> h
        (
            const scalarField& T,
            const labelList& cells
        ) const;


        //- Enthalpy for patch [J/kg]
        virtual tmp<scalarField> h
        (
            const scalarField& T,
            const label patchi
        ) const;


        //- Density for patch [J/kg]
        virtual tmp<scalarField> rho
        (
            const scalarField& p,
            const scalarField& h,
            const label patchi
        ) const;

        //- Heat capacity at constant pressure for patch [J/kg/K]
        // dummy function needed for BC
        virtual tmp<scalarField> Cp
        (
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        virtual tmp<volScalarField> Cp() const;

        //- Heat capacity at constant volume for patch [J/kg/K]
        virtual tmp<scalarField> Cv
        (
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        virtual tmp<volScalarField> Cv() const;

        //- Gradient drhodh @ constant pressure
        virtual const volScalarField& drhodh() const
        {
            return drhodh_;
        }

        //- Density [kg/m^3] - uses current value of pressure
        virtual tmp<volScalarField> rho() const;

        //- Compressibility [s^2/m^2]
        virtual const volScalarField& psi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IAPWSThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
