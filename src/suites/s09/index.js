// Suite 9 — Risers. Modules are listed in the order of the specified analyses.
import staticConfiguration from './static-configuration.js';
import coupledDynamic from './coupled-dynamic.js';
import tensionMoment from './tension-moment.js';
import vesselMotion from './vessel-motion.js';
import viv from './viv.js';
import waveFatigue from './wave-fatigue.js';
import touchdownZone from './touchdown-zone.js';
import bucklingCollapse from './buckling-collapse.js';
import interference from './interference.js';
import hangoffLoads from './hangoff-loads.js';
import installation from './installation.js';
import internalFlow from './internal-flow.js';
import extremeDisconnect from './extreme-disconnect.js';
import fatigueIntegrity from './fatigue-integrity.js';
import localStrength from './local-strength.js';

export default {
  n: 9,
  modules: [staticConfiguration, coupledDynamic, tensionMoment, vesselMotion, viv, waveFatigue, touchdownZone, bucklingCollapse, interference, hangoffLoads,
    installation, internalFlow, extremeDisconnect, fatigueIntegrity, localStrength],
};
