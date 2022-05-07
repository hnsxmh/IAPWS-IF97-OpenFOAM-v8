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

class IAPWSThermo
:
    public fluidThermo
{
    // Private data

    //- DensityField
        volScalarField he_;

        volScalarField rho_;

        volScalarField drhodh_;

        volScalarField psi_;

        //- Laminar dynamic viscosity [kg/m/s]
        volScalarField mu_;

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

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);

        //- Disallow default bitwise copy construction
        IAPWSThermo(const IAPWSThermo&) = delete;


    //- Destructor
    virtual ~IAPWSThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const IAPWSThermo&
        ) = delete;

    // Fields derived from thermodynamic state variables

        //- Enthalpy for cell-set [J/kg]
        virtual tmp<scalarField> he
        (
            const scalarField& T,
            const labelList& cells
        ) const;


        //- Enthalpy for patch [J/kg]
        virtual tmp<scalarField> he
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

        // //- Return non-const access to the local density field [kg/m^3]
        // virtual volScalarField& rho()
        // {
        //     return rho_;
        // }

        //- Compressibility [s^2/m^2]
        virtual const volScalarField& psi() const;

        //- Enthalpy/Internal energy [J/kg]
        //  Non-const access allowed for transport equations
        virtual volScalarField& he()
        {
            return he_;
        }

        //- Enthalpy/Internal energy [J/kg]
        virtual const volScalarField& he() const
        {
            return he_;
        }

        virtual tmp<volScalarField> mu() const
        {
            return mu_;
        }

        virtual tmp<scalarField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }
        
        virtual tmp<scalarField> rho(const label patchi) const
        {
            return rho_.boundaryField()[patchi];
        }

//useless func from base class
        virtual tmp<volScalarField> kappaEff
        (
            const volScalarField& alphat
        ) const
        {
            return volScalarField::New
            (
                "kappaEff",
                Cp()*(this->alpha_ + alphat)
            );
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return
                Cp
                (
                    this->T_.boundaryField()[patchi],
                    patchi
                )
            *(this->alpha_.boundaryField()[patchi] + alphat);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const
        {
            return volScalarField::New
            (
                "alphaEff",
                this->alpha_ + alphat
            );
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return this->alpha_.boundaryField()[patchi] + alphat;
        }

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return word("IAPWSThermo");
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return false;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return false;
        }

        //- Enthalpy/Internal energy
        //  for given pressure and temperature [J/kg]
        virtual tmp<volScalarField> he
        (
            const volScalarField& p,
            const volScalarField& T
        ) const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Sensible enthalpy [J/kg]
        virtual tmp<volScalarField> hs() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Sensible enthalpy
        //  for given pressure and temperature [J/kg]
        virtual tmp<volScalarField> hs
        (
            const volScalarField& p,
            const volScalarField& T
        ) const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Sensible enthalpy for cell-set [J/kg]
        virtual tmp<scalarField> hs
        (
            const scalarField& T,
            const labelList& cells
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Sensible enthalpy for patch [J/kg]
        virtual tmp<scalarField> hs
        (
            const scalarField& T,
            const label patchi
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Absolute enthalpy [J/kg]
        virtual tmp<volScalarField> ha() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Absolute enthalpy
        //  for given pressure and temperature [J/kg]
        virtual tmp<volScalarField> ha
        (
            const volScalarField& p,
            const volScalarField& T
        ) const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Absolute enthalpy for cell-set [J/kg]
        virtual tmp<scalarField> ha
        (
            const scalarField& T,
            const labelList& cells
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Absolute enthalpy for patch [J/kg]
        virtual tmp<scalarField> ha
        (
            const scalarField& T,
            const label patchi
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Enthalpy of formation [J/kg]
        virtual tmp<volScalarField> hc() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Temperature from enthalpy/internal energy for cell-set
        virtual tmp<scalarField> THE
        (
            const scalarField& h,
            const scalarField& T0,      // starting temperature
            const labelList& cells
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Temperature from enthalpy/internal energy for patch
        virtual tmp<scalarField> THE
        (
            const scalarField& h,
            const scalarField& T0,      // starting temperature
            const label patchi
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Heat capacity at constant pressure/volume [J/kg/K]
        virtual tmp<volScalarField> Cpv() const
        {
            return this->Cp();
        }

        //- Heat capacity at constant pressure/volume for patch [J/kg/K]
        virtual tmp<scalarField> Cpv
        (
            const scalarField& T,
            const label patchi
        ) const
        {
            return Cp(T,patchi);
        }

        //- Heat capacity ratio []
        virtual tmp<volScalarField> CpByCpv() const
        {
            const fvMesh& mesh = this->T_.mesh();

            tmp<volScalarField> tCpByCpv
            (
                new volScalarField
                (
                    IOobject
                    (
                        "CpByCpv",
                        mesh.time().timeName(),
                        this->T_.db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("", dimensionSet(0, 0, 0, 0, 0), 1.0)
                )
            );

            return tCpByCpv;
        }

        //- Heat capacity ratio for patch []
        virtual tmp<scalarField> CpByCpv
        (
            const scalarField& T,
            const label patchi
        ) const
        {
            tmp<scalarField> tPsi
            (
                new scalarField(this->T_.boundaryField()[patchi].size())
            );
            scalarField& psi = tPsi.ref();

            forAll(this->T_.boundaryField()[patchi], facei)
            {
                psi[facei] = 1.0;
            }

            return tPsi;
        }

        //- Thermal diffusivity for temperature of mixture [W/m/K]
        virtual tmp<volScalarField> kappa() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Thermal diffusivity for temperature of mixture
        //  for patch [W/m/K]
        virtual tmp<scalarField> kappa
        (
            const label patchi
        ) const
        {
            return tmp<scalarField>(nullptr);
        }

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const
        {
            return volScalarField::New
            (
                "alphahe",
                this->alpha_
            );
        }

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const
        {
            tmp<scalarField> tPsi
            (
                new scalarField(this->T_.boundaryField()[patchi].size())
            );
            scalarField& psi = tPsi.ref();

            forAll(this->T_.boundaryField()[patchi], facei)
            {
                psi[facei] = this->alpha_.boundaryField()[patchi][facei];
            }

            return tPsi;
        }

        virtual tmp<volScalarField> W() const
        {
            return tmp<volScalarField>(nullptr);
        }

        virtual tmp<scalarField> W(const label patchi) const
        {
            return tmp<scalarField>(nullptr);
        }

        virtual tmp<volScalarField> gamma() const
        {
            return tmp<volScalarField>(nullptr);
        }

        virtual tmp<scalarField> gamma(const scalarField&, label) const
        {
            return tmp<scalarField>(nullptr);
        }

        virtual void correctRho(const volScalarField& deltaRho)
        {}
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
