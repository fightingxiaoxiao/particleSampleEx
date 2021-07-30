#ifndef SampleParticle_H
#define SampleParticle_H

#include "fvMesh.H"

namespace Foam
{
    class SampleParticle
    {
    public:
        label id;

        scalar d;

        scalar rho;

        label nParticle;

        vector position;

        vector U;

        scalar mass();

        SampleParticle();

        SampleParticle(scalar _d,
                       scalar _rho,
                       label _nParticle,
                       vector _position,
                       vector _U);

        ~SampleParticle();
    };
}

#endif
