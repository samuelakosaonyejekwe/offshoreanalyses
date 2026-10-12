import globalResponse from './global-response.js';
import crossSection from './cross-section.js';
import tensionBending from './tension-bending.js';
import torsion from './torsion.js';
import burst from './burst.js';
import collapse from './collapse.js';
import wetCollapse from './wet-collapse.js';
import armourFatigue from './armour-fatigue.js';
import birdcaging from './birdcaging.js';
import mbr from './mbr.js';
import interlayerContact from './interlayer-contact.js';
import thermalPressure from './thermal-pressure.js';
import annulus from './annulus.js';
import endFitting from './end-fitting.js';
import installation from './installation.js';
import serviceLife from './service-life.js';

export default {
  n: 14,
  modules: [globalResponse, crossSection, tensionBending, torsion, burst, collapse, wetCollapse, armourFatigue, birdcaging, mbr, interlayerContact, thermalPressure, annulus, endFitting, installation, serviceLife],
};
