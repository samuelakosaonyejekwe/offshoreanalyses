// Suite 11 — Floating Offshore Structures.
import hydrostatics from './hydrostatics.js';
import intactStability from './intact-stability.js';
import waveLoads from './wave-loads.js';
import motions from './motions-6dof.js';
import rao from './rao.js';
import frequencyDomain from './frequency-domain.js';
import timeDomain from './time-domain.js';
import globalStrength from './global-strength.js';
import localStrength from './local-strength.js';
import fatigue from './fatigue.js';
import coupled from './coupled-mooring-riser.js';
import survival from './survival.js';
import damageStability from './damage-stability.js';
import progressiveFlooding from './progressive-flooding.js';
import slamming from './slamming-green-water.js';
import vim from './vim.js';
import ballast from './ballast.js';
import dynamicPositioning from './dynamic-positioning.js';
import accidentalDamage from './accidental-damage.js';
import integrityLife from './integrity-life.js';

export default {
  n: 11,
  modules: [hydrostatics, intactStability, waveLoads, motions, rao, frequencyDomain, timeDomain, globalStrength, localStrength, fatigue, coupled,
    survival, damageStability, progressiveFlooding, slamming, vim, ballast, dynamicPositioning, accidentalDamage, integrityLife],
};
