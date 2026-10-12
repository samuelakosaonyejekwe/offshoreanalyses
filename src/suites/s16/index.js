import coupledDynamics from './coupled-dynamics.js';
import wecHydrodynamics from './wec-hydrodynamics.js';
import tidalTurbine from './tidal-turbine.js';
import rotorAerodynamics from './rotor-aerodynamics.js';
import globalStrength from './global-strength.js';
import fatigue from './fatigue.js';
import floaterStability from './floater-stability.js';
import mooring from './mooring.js';
import foundationSoil from './foundation-soil.js';
import extremeLoads from './extreme-loads.js';
import naturalFrequency from './natural-frequency.js';
import ptoDynamics from './pto-dynamics.js';
import generatorConverter from './generator-converter.js';
import energyProduction from './energy-production.js';
import floatingSolar from './floating-solar.js';
import wakeArray from './wake-array.js';
import arrayCable from './array-cable.js';
import installationMaintenance from './installation-maintenance.js';
import scour from './scour.js';
import reliabilityLifetime from './reliability-lifetime.js';

// Ordered as the analyses are listed in the suite specification; the generator and converter chain follows the
// power-take-off dynamics, and floating solar follows the energy-yield module it shares its analysis with.
export default {
  n: 16,
  modules: [coupledDynamics, wecHydrodynamics, tidalTurbine, rotorAerodynamics, globalStrength, fatigue, floaterStability, mooring, foundationSoil, extremeLoads, naturalFrequency, ptoDynamics,
    generatorConverter, energyProduction, floatingSolar, wakeArray, arrayCable, installationMaintenance, scour, reliabilityLifetime],
};
