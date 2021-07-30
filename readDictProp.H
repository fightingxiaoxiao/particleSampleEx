IOdictionary propsDict(
    IOobject(
        "particleSampleProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED));

const word sampleDirection(propsDict.get<word>("sampleDirection"));

const scalar samplePosition(propsDict.get<scalar>("samplePosition"));

const scalar limitMoveDistanceInOneSample(propsDict.get<scalar>("limitMoveDistanceInOneSample"));

const scalar startHeight(propsDict.get<scalar>("startHeight"));

const scalar deltaH(propsDict.get<scalar>("deltaH"));

label dynamicParticleListLength(propsDict.get<label>("dynamicParticleListLength"));

label directionIndex;
if (sampleDirection == "x")
{
    directionIndex = 0;
}
else if (sampleDirection == "y")
{
    directionIndex = 1;
}
else if (sampleDirection == "z")
{
    directionIndex = 2;
}
else
{
    Info << "Error direction name." << endl;
    return 1;
}
