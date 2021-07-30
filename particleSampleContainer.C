#include "particleSampleContainer.H"

namespace Foam
{

    void particleSampleContainer::classifyDiameterAlongHeight(scalar startHeight,
                                                              scalar deltaH)
    {
        for (auto &p : particleStorage)
        {
            if (p.second.position[2] - startHeight < 0)
                continue;
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);
            diameterList[hIndex].push_back(p.second.d);
        }
    }

    void particleSampleContainer::classifyVelocityAlongHeight(scalar startHeight,
                                                              scalar deltaH)
    {
        for (auto &p : particleStorage)
        {
            if (p.second.position[2] - startHeight < 0)
                continue;
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);
            velocityList[hIndex].push_back(mag(p.second.U));
        }
    }

    void particleSampleContainer::classifyFlowRateAlongHeight(scalar startHeight,
                                                              scalar deltaH,
                                                              scalar directionIndex,
                                                              scalar samplePosition,
                                                              scalar limitMoveDistanceInOneSample)
    {
        for (auto &p : particleStorage)
        {
            if (p.second.position[2] - startHeight < 0)
                continue;
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);

            if (massRateList.find(hIndex) == massRateList.end())
            {
                massRateList[hIndex] = 0.;
            }

            auto pOld = _particleStorage[p.first];
            auto _position = pOld.position;
            auto position = p.second.position;
            if (mag(position - _position) > limitMoveDistanceInOneSample)
            {
                continue;
            }

            if (_position[directionIndex] < samplePosition && position[directionIndex] >= samplePosition)
            {
                massRateList[hIndex] += p.second.mass();
            }
            else if (_position[directionIndex] > samplePosition && position[directionIndex] <= samplePosition)
            {
                massRateList[hIndex] -= p.second.mass();
            }
        }
    }

    void particleSampleContainer::classifyFlowRateAlongHeight(scalar startHeight,
                                                              scalar deltaH)
    {
        for (auto &p : particleStorage)
        {
            if (p.second.position[2] - startHeight < 0)
                continue;
            label hIndex = static_cast<label>((p.second.position[2] - startHeight) / deltaH);

            if (massRateList.find(hIndex) == massRateList.end())
            {
                massRateList[hIndex] = 0.;
            }

            massRateList[hIndex] += p.second.mass() * mag(p.second.U);
        }
    }

    word particleSampleContainer::writeDiameterInfo(scalar startHeight,
                                                    scalar deltaH)
    {
        std::stringstream diameterInfo("");
        diameterInfo << "height, diameter0, diameter1..." << std::endl;
        for (auto &d : diameterList)
        {
            diameterInfo << startHeight + d.first * deltaH;

            for (auto &data : d.second)
            {
                diameterInfo << "," << data;
            }

            diameterInfo << std::endl;
        }
        return diameterInfo.str();
    }

    word particleSampleContainer::writeVelocityInfo(scalar startHeight,
                                                    scalar deltaH)
    {
        std::stringstream velocityInfo("");
        velocityInfo << "height, velmag0, velmag1..." << std::endl;
        for (auto &v : velocityList)
        {
            velocityInfo << startHeight + v.first * deltaH;

            for (auto &data : v.second)
            {
                velocityInfo << "," << data;
            }

            velocityInfo << std::endl;
        }

        return velocityInfo.str();
    }

    word particleSampleContainer::writeFlowRateInfo(scalar startHeight,
                                                    scalar deltaH)
    {
        std::stringstream writeFlowRateInfo("");
        for (auto &m : massRateList)
        {
            writeFlowRateInfo << startHeight + m.first * deltaH << ", " << m.second << std::endl;
        }

        return writeFlowRateInfo.str();
    }

    scalar particleSampleContainer::writeTotalFlowRate(scalar startHeight)
    {
        scalar totalFlowRate = 0.;
        for (auto &m : massRateList)
        {
            totalFlowRate += m.second;
        }

        return totalFlowRate;
    }

}
