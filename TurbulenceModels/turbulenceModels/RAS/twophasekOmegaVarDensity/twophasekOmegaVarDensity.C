/*---------------------------------------------------------------------------*\
Copyright (C) 2015 Cyrille Bonamy, Julien Chauchat, Tian-Jian Hsu
                   and contributors

License
    This file is part of SedFOAM.

    SedFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SedFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with SedFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twophasekOmegaVarDensity.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void twophasekOmegaVarDensity<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = k_/
    (
        max(omega_, Clim_*sqrt((2.0*(magSqr(symm(fvc::grad(this->U_)))))/Cmu_))
    );

    this->nut_.min(nutMax_);
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
twophasekOmegaVarDensity<BasicTurbulenceModel>::twophasekOmegaVarDensity
(
    const alphaField& beta,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& betaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        beta,
        rho,
        U,
        betaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    popeCorrection_
    (
        Switch::getOrAddToDict
        (
            "popeCorrection",
            this->coeffDict_,
            true
        )
    ),
    writeTke_
    (
        Switch::getOrAddToDict
        (
            "writeTke",
            this->coeffDict_,
            false
        )
    ),
    C3om_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C3om",
            this->coeffDict_,
            0.35
        )
    ),
    C4om_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C4om",
            this->coeffDict_,
            1.0
        )
    ),
    KE2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE2",
            this->coeffDict_,
            1.0
        )
    ),
    KE4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "KE4",
            this->coeffDict_,
            1.0
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    betaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaOmega",
            this->coeffDict_,
            0.072
        )
    ),
    nutMax_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "nutMax",
            this->coeffDict_,
            1e-1
        )
    ),
    Clim_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Clim",
            this->coeffDict_,
            0.0
        )
    ),
    sigmad_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmad",
            this->coeffDict_,
            0.0
        )
    ),
    alphaKOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaKOmega",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.52
        )
    ),
    alphaOmegaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmegaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    /*sigmat_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmat",
            this->coeffDict_,
            0.85
        )
    ),*/
    tmfexp_(U.db().lookupObject<volScalarField> ("tmfexp")),
    ESD3_(U.db().lookupObject<volScalarField> ("ESD3")),
    ESD4_(U.db().lookupObject<volScalarField> ("ESD4")),
    ESD5_(U.db().lookupObject<volScalarField> ("ESD5")),
    ESD_(U.db().lookupObject<volScalarField> ("ESD")),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool twophasekOmegaVarDensity<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        betaOmega_.readIfPresent(this->coeffDict());
        alphaOmegaOmega_.readIfPresent(this->coeffDict());
        alphaKOmega_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void twophasekOmegaVarDensity<BasicTurbulenceModel>::correct()
{
    if (not this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& beta = this->alpha_;
    const alphaField alpha = 1 - beta;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& betaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const surfaceScalarField& phi = this->phi_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    volSymmTensorField Sij(symm(fvc::grad(U)));

    volScalarField G
    (
        this->GName(),
        nut*2*magSqr(Sij)
    );

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volTensorField Omij(-skew(fvc::grad(U)));
    volVectorField Gradk(fvc::grad(k_));
    volVectorField Gradomega(fvc::grad(omega_));
    volScalarField alphadCheck_(Gradk & Gradomega);

    const volScalarField CDkOmega
    (
    sigmad_*pos(0.15-alpha)*pos(alphadCheck_)*(alphadCheck_)/omega_
    );


    volScalarField XsiOmega
    (
        IOobject
        (
            IOobject::groupName("XsiOmega", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    if (popeCorrection_)
    {
        XsiOmega =
        (
            mag((Omij & Omij & Sij)/(pow((Cmu_*omega_), 3)))
        );
    }

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(beta, rho, omega_)
      + fvm::div(betaRhoPhi, omega_)
      - fvm::Sp(fvc::div(phi)*beta*rho, omega_)
      - fvm::laplacian(beta*rho*DomegaEff(), omega_)
      ==
      - fvm::SuSp (-beta()*rho()*alphaOmega_*G/k_, omega_)
      - fvm::Sp(beta()*rho()*ESD_, omega_)
      - fvm::Sp
        (
            beta()*rho()*betaOmega_*
            (
                (scalar(1.0)+scalar(85.0)*XsiOmega())
                /(scalar(1.0)+scalar(100.0)*XsiOmega())
            )*omega_(),
            omega_
        )
      + beta()*rho()*CDkOmega
      + ESD2()*fvm::Sp(beta()*rho()*C3om_*KE2_, omega_)
      + fvm::Sp((beta()*rho()*C4om_*KE4_*ESD5_*nut/k_), omega_)
    );
    if (writeTke_)
    {
        #include "writeTKEBudget_kOmega.H"
    }
    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(beta, rho, k_)
      + fvm::div(betaRhoPhi, k_)
      - fvm::Sp(fvc::div(phi)*beta*rho, k_)
      - fvm::laplacian(beta*rho*DkEff(), k_)
      ==
      - fvm::SuSp(-beta()*rho()*G/k_, k_)
      + fvm::Sp(-beta()*rho()*Cmu_*omega_, k_)
      + fvm::Sp(beta()*rho()*ESD_, k_)
      + fvm::Sp(beta()*rho()*KE4_*ESD4_*nut/k_, k_)
      + ESD2()*fvm::Sp(beta()*rho()*KE2_, k_)
      //- fvm::SuSp((g & fvc::grad(rho))/(omega_*sigmat), k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
