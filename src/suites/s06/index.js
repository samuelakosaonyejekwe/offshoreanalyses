// Suite 6 — Dolphins, mooring dolphins and berthing dolphins.
import berthingEnergy from './berthing-energy.js';
import fenderPerformance from './fender-performance.js';
import berthingDynamics from './berthing-dynamics.js';
import mooringLines from './mooring-lines.js';
import environmentalLoading from './environmental-loading.js';
import globalStrength from './global-strength.js';
import pileCapacity from './pile-capacity.js';
import soilPileInteraction from './soil-pile-interaction.js';
import flexibleDolphin from './flexible-dolphin.js';
import bucklingStability from './buckling-stability.js';
import fatigueImpact from './fatigue-impact.js';
import seismicResponse from './seismic-response.js';
import integrityDamage from './integrity-damage.js';

export default {
  n: 6,
  modules: [berthingEnergy, fenderPerformance, berthingDynamics, mooringLines, environmentalLoading, globalStrength, pileCapacity, soilPileInteraction, flexibleDolphin, bucklingStability, fatigueImpact, seismicResponse, integrityDamage],
};
