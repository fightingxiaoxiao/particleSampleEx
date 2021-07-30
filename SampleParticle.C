#include "SampleParticle.H"

namespace Foam
{

    scalar SampleParticle::mass()
    {
        return nParticle * 4 / 3 * M_PI * Foam::pow(d / 2., 3.) * rho;
    }

    SampleParticle::SampleParticle(scalar _d,
                                   scalar _rho,
                                   label _nParticle,
                                   vector _position,
                                   vector _U)
    {
        d = _d;
        rho = _rho;
        nParticle = _nParticle;
        position = _position;
        U = _U;
    }

    SampleParticle::SampleParticle() {}

    SampleParticle::~SampleParticle() {}
}
