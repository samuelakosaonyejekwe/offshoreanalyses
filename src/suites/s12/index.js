import staticConfig from './static-config.js';
import coupledDynamics from './coupled-dynamics.js';
import tensionBending from './tension-bending.js';
import mbr from './mbr.js';
import crossSection from './cross-section.js';
import fatigue from './fatigue.js';
import seabedContact from './seabed-contact.js';
import installation from './installation.js';
import vesselMotion from './vessel-motion.js';
import viv from './viv.js';
import crushCollapse from './crush-collapse.js';
import thermal from './thermal.js';
import tubePressure from './tube-pressure.js';
import interference from './interference.js';
import serviceLife from './service-life.js';

// Ordered as the analyses are listed in the suite specification.
export default {
  n: 12,
  modules: [staticConfig, coupledDynamics, tensionBending, mbr, crossSection, fatigue, seabedContact, installation, vesselMotion, viv, crushCollapse, thermal, tubePressure, interference, serviceLife],
};
