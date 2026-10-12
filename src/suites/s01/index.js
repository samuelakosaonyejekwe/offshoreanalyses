// Suite 1 — Fixed offshore platforms and jackets.
import waveLoading from './wave-loading.js';
import inplace from './inplace.js';
import staticStrength from './static-strength.js';
import ulsSls from './uls-sls.js';
import globalLocal from './global-local.js';
import joints from './joints.js';
import stability from './stability.js';
import nonlinear from './nonlinear.js';
import pushover from './pushover.js';
import progressiveCollapse from './progressive-collapse.js';
import dynamics from './dynamics.js';
import seismic from './seismic.js';
import fatigue from './fatigue.js';
import pileFoundation from './pile-foundation.js';
import collision from './collision.js';
import loadoutTransport from './loadout-transport.js';
import lifting from './lifting.js';
import reassessment from './reassessment.js';

export default {
  n: 1,
  modules: [
    waveLoading, inplace, staticStrength, ulsSls, globalLocal, joints, stability, nonlinear, pushover, progressiveCollapse,
    dynamics, seismic, fatigue, pileFoundation, collision, loadoutTransport, lifting, reassessment,
  ],
};
