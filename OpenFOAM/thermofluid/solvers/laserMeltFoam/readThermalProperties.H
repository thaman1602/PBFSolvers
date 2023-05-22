IOdictionary thermalProperties
(
    IOobject
    (
    	"thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
     )
);

dictionary metal = thermalProperties.subDict("metalProperties");
dictionary gas = thermalProperties.subDict("gasProperties");

Info<< "\nReading thermalProperties\n" << endl;
//const dimensionedScalar betaM (metal.lookup("betaM"));
const dimensionedScalar Tl ("Tl", dimensionSet(0,0,0,1,0,0,0),metal.lookup("Tl"));
const dimensionedScalar Ts ("Ts", dimensionSet(0,0,0,1,0,0,0),metal.lookup("Ts"));
const dimensionedScalar Hf ("Hf", dimensionSet(0,2,-2,0,0,0,0),metal.lookup("Hf"));
const dimensionedScalar C ("C", dimensionSet(1,-3,-1,0,0,0,0),metal.lookup("C"));
//const dimensionedScalar betaM ("betaM", dimensionSet(0,0,0,-1,0,0,0),metal.lookup("betaM"));
const dimensionedVector kSCoeff ("kSCoeff",dimensionSet(1,1,-3,-1,0,0,0),metal.lookup("kSCoeff")); //Coefficients of the solid thermal conductivity polynomial
const dimensionedVector kLCoeff ("kLCoeff",dimensionSet(1,1,-3,-1,0,0,0),metal.lookup("kLCoeff")); //Coefficients of the solid thermal conductivity polynomial
const dimensionedVector cpSCoeff ("cpSCoeff",dimensionSet(0,2,-2,-1,0,0,0),metal.lookup("cpSCoeff"));
const dimensionedScalar cpL ("cpL", dimensionSet(0,2,-2,-1,0,0,0),metal.lookup("cpL"));

const dimensionedScalar kG ("kG",dimensionSet(1,1,-3,-1,0,0,0),gas.lookup("kG")); //Constant gas thermal conductivity
const dimensionedScalar cpG ("cpG", dimensionSet(0,2,-2,-1,0,0,0),gas.lookup("cpG"));

const dimensionedVector rhoLCoeff ("rhoLCoeff",dimensionSet(1,-3,0,0,0,0,0),metal.lookup("rhoLCoeff"));
//Vapour Pressure parameters
Info<<"Reading vapour pressure parameters \n" <<endl;

const dimensionedScalar p0 ("p0", dimensionSet(1,-1,-2,0,0,0,0),metal.lookup("p0")); // Atmospheric pressure
const dimensionedScalar TVap ("TVap", dimensionSet(0,0,0,1,0,0,0),metal.lookup("TVap")); // Vaporisation temperature
const dimensionedScalar Hv ("Hv", dimensionSet(0,2,-2,0,0,0,0),metal.lookup("Hv")); // Enthalpy change due to vaporisation
const dimensionedScalar R ("R", dimensionSet(1,2,-2,-1,-1,0,0),metal.lookup("R")); // Universal gas constant
const dimensionedScalar Mm ("Mm", dimensionSet(1,0,0,0,-1,0,0),metal.lookup("Mm")); //Molar mass of metal

const dimensionedScalar TRef ("TRef", dimensionSet(0,0,0,1,0,0,0),thermalProperties.lookup("TRef"));
const dimensionedScalar hC ("hC", dimensionSet(1,0,-3,-1,0,0,0),thermalProperties.lookup("hC"));


dimensionedScalar TSmooth("TSmooth",dimensionSet(0, 0, 0, 1, 0),50.0);
