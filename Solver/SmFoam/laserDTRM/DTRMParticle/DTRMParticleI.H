/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::DTRMParticle::trackingData::trackingData
(
    Cloud<DTRMParticle>& spc,
    const interpolationCellPoint<vector>& nHatInterp,
    const labelField& relfectedCell,
    //const UPtrList<reflectionModel>& reflection,
    volScalarField& Q
)
:
    particle::trackingData(spc),
    nHatInterp_(nHatInterp),
    relfectedCells_(relfectedCell),
    //reflection_(reflection),
    Q_(Q)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::DTRMParticle::trackingData::nHatInterp() const
{
    return nHatInterp_;
}

inline const Foam::labelField&
Foam::DTRMParticle::trackingData::relfectedCells() const
{
    return relfectedCells_;
}

/*
inline const Foam::UPtrList<Foam::radiation::reflectionModel>&
Foam::DTRMParticle::trackingData::reflection() const
{
    return reflection_;
}

*/
inline Foam::scalar& Foam::DTRMParticle::trackingData::Q(label celli)
{
    return Q_[celli];
}


inline const Foam::point& Foam::DTRMParticle::p0() const
{
    return p0_;
}


inline const Foam::point& Foam::DTRMParticle::p1() const
{
    return p1_;
}


inline Foam::scalar Foam::DTRMParticle::I0() const
{
    return I0_;
}


inline Foam::scalar Foam::DTRMParticle::I() const
{
    return I_;
}


inline Foam::scalar Foam::DTRMParticle::dA() const
{
    return dA_;
}


inline Foam::scalar& Foam::DTRMParticle::dA()
{
    return dA_;
}


inline Foam::point& Foam::DTRMParticle::p1()
{
    return p1_;
}


inline Foam::scalar& Foam::DTRMParticle::I0()
{
    return I0_;
}


inline Foam::scalar& Foam::DTRMParticle::I()
{
    return I_;
}


// ************************************************************************* //
