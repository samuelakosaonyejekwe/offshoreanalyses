import staticConfig from './static-config.js';
import dynamicResponse from './dynamic-response.js';
import tensionBendingTorsion from './tension-bending-torsion.js';
import fatigue from './fatigue.js';
import touchdownSeabed from './touchdown-seabed.js';
import layingInstallation from './laying-installation.js';
import burialProtection from './burial-protection.js';
import thermalRating from './thermal-rating.js';
import electricalLosses from './electrical-losses.js';
import opticalLink from './optical-link.js';
import emfInduced from './emf-induced.js';
import crushingPressure from './crushing-pressure.js';
import freeSpanViv from './free-span-viv.js';
import crossingInterference from './crossing-interference.js';
import failureReliability from './failure-reliability.js';
import extremeStorm from './extreme-storm.js';

// Ordered as the analyses are listed in the suite specification; the optical link follows the electrical loading it belongs with.
export default {
  n: 13,
  modules: [staticConfig, dynamicResponse, tensionBendingTorsion, fatigue, touchdownSeabed, layingInstallation, burialProtection, thermalRating, electricalLosses, opticalLink,
    emfInduced, crushingPressure, freeSpanViv, crossingInterference, failureReliability, extremeStorm],
};
