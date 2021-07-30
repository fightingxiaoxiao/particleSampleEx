#ifndef particleSampleContainer_H
#define particleSampleContainer_H

#include <vector>
#include <map>
#include "SampleParticle.H"

namespace Foam
{
    class particleSampleContainer
    {
    private:
    public:
        std::map<label, SampleParticle, std::less<label>> particleStorage;

        std::map<label, SampleParticle, std::less<label>> _particleStorage;

        std::map<label, std::vector<scalar>, std::less<label>> diameterList;

        std::map<label, std::vector<scalar>, std::less<label>> velocityList;

        std::map<label, scalar, std::less<label>> massRateList;

        void classifyDiameterAlongHeight(scalar startHeight, scalar deltaH);

        void classifyVelocityAlongHeight(scalar startHeight, scalar deltaH);

        void classifyFlowRateAlongHeight(scalar startHeight,
                                         scalar deltaH,
                                         scalar directionIndex,
                                         scalar samplePosition,
                                         scalar limitMoveDistanceInOneSample);

        void classifyFlowRateAlongHeight(scalar startHeight,
                                         scalar deltaH);

        word writeDiameterInfo(scalar startHeight,
                               scalar deltaH);

        word writeVelocityInfo(scalar startHeight,
                               scalar deltaH);

        word writeFlowRateInfo(scalar startHeight,
                               scalar deltaH);

        scalar writeTotalFlowRate(scalar startHeight);
    };
}
#endif
