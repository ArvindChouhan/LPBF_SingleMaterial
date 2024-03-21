/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd
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

#include "DTRMParticle.H"
#include "constants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    const vector& position,
    const vector& targetPosition,
    const scalar I,
    const label cellI,
    const scalar dA,
    const label transmissiveId
)
:
    particle(mesh, position, cellI),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    transmissiveId_(transmissiveId)
{}


Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const vector& position,
    const vector& targetPosition,
    const scalar I,
    const scalar dA,
    const label transmissiveId
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    transmissiveId_(transmissiveId)
{}


Foam::DTRMParticle::DTRMParticle(const DTRMParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    I0_(p.I0_),
    I_(p.I_),
    dA_(p.dA_),
    transmissiveId_(p.transmissiveId_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::DTRMParticle::move
(
    Cloud<DTRMParticle>& spc,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        //Cache old data of particle to use for reflected particle
        const point pos0 = position();
        const label cell1 = cell();    //cell 1

        scalar f = 1 - stepFraction();
        const vector s = p1() - p0() - deviationFromMeshCentre();
        trackToAndHitFace(f*s, f, spc, td);

        const point p1 = position();
        vector dsv = p1 - pos0;
        scalar ds = mag(dsv);

        //const label cell1 = cell();	//Delete Comment

        //NOTE:
        // Under the new barocentric tracking alghorithm the newly
        // inserted particles are tracked to the nearest cell centre first,
        // then, given the direction, to a face. In both occasions the first call
        // to trackToAndHitFace returns ds = 0. In this case we do an extra
        // call to trackToAndHitFace to start the tracking.
        // This is a temporary fix until the tracking can handle it.
        if (ds == 0)
        {
            trackToAndHitFace(f*s, f, spc, td);
            dsv = p1 - position();
            ds = mag(dsv);
        }

        scalar myI0_ = 17.897e9;

        label reflectedZoneId = td.relfectedCells()[cell1];

        if
        (
            (reflectedZoneId > -1)
         && (
                (transmissiveId_ == -1)
             || (transmissiveId_ != reflectedZoneId)
            )
        )
        {
            scalar rho(0);

            // Create a new reflected particle when the particles is not
            // transmissive and larger than an absolute I
            if (I_ > 0.01*myI0_ && ds > 0)  //I0_
            {
                vector pDir = dsv/ds;

                cellPointWeight cpw(mesh(), pos0, cell1, face()); // position()
                vector nHat = td.nHatInterp().interpolate(cpw);

                nHat /= (mag(nHat) + ROOTSMALL);
                scalar cosTheta(-pDir & nHat);
		scalar epsilon_ = 0.0625;             //0.0625
                // Only new incoming rays
                if (cosTheta > SMALL)
                {
                    vector newDir = pDir + 2.0*(-pDir & nHat) * nHat;

                    // reflectivity
                    rho =
                        min
                        (
                            max
                            (
        0.5
      * (
            (1 + sqr(1 - epsilon_*cosTheta))/(1 + sqr(1 + epsilon_*cosTheta))
        +
            (sqr(epsilon_) - 2*epsilon_*cosTheta + 2*sqr(cosTheta))
          /
            (sqr(epsilon_) + 2*epsilon_*cosTheta + 2*sqr(cosTheta))
        ), 0.0
                            )
                            , 0.98
                        );

                    scalar delaM = cbrt(mesh().cellVolumes()[cell1]);

                    const point insertP(pos0 - pDir*0.1*delaM);  //  position()
                    label cellI = mesh().findCell(insertP);

                    if (cellI > -1)
                    {
                        DTRMParticle* pPtr = new DTRMParticle
                        (
                            mesh(),
                            insertP,
                            insertP + newDir*mesh().bounds().mag(),
                            I_*rho,
                            cellI,
                            dA_,
                            -1
                        );

                        // Add to cloud
                        spc.addParticle(pPtr);
                    }

            transmissiveId_ = reflectedZoneId;

            const scalar Itran = I_*(1.0 - rho);

            td.Q(cell1) += (Itran)*dA_;

	    break;
                }
            }
        }
        else
        {
            if ((I_ <= 0.01*myI0_))  //I0_
            {
                stepFraction() = 1.0;
                break;
            }
        }

    }

    return td.keepParticle;
}


void Foam::DTRMParticle::hitProcessorPatch
(
    Cloud<DTRMParticle>&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::DTRMParticle::hitWallPatch
(
    Cloud<DTRMParticle>&,
    trackingData& td
)
{
    td.keepParticle = false;
}


bool Foam::DTRMParticle::hitPatch
(
    Cloud<DTRMParticle>&,
    trackingData& td
)
{
    return false;
}


// ************************************************************************* //
