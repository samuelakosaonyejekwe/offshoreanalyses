// Suite 3 — Offshore wind turbine foundations: monopiles, jackets, tripods, gravity bases and suction buckets.
import coupledSimulation from './coupled-simulation.js';
import structuralStrength from './structural-strength.js';
import windWaveCurrentLoading from './wind-wave-current-loading.js';
import dynamicResponse from './dynamic-response.js';
import fatigueLife from './fatigue-life.js';
import monopileCapacity from './monopile-capacity.js';
import pileSoilInteraction from './pile-soil-interaction.js';
import gravityBaseStability from './gravity-base-stability.js';
import bucklingUltimate from './buckling-ultimate.js';
import scour from './scour.js';
import seismicLiquefaction from './seismic-liquefaction.js';
import naturalFrequency from './natural-frequency.js';
import extremeLoading from './extreme-loading.js';
import pileDriving from './pile-driving.js';
import suctionBucketInstallation from './suction-bucket-installation.js';
import vesselImpact from './vessel-impact.js';
import reliabilityLifeExtension from './reliability-life-extension.js';

export default {
  n: 3,
  modules: [
    coupledSimulation, structuralStrength, windWaveCurrentLoading, dynamicResponse, fatigueLife, monopileCapacity, pileSoilInteraction, gravityBaseStability, bucklingUltimate,
    scour, seismicLiquefaction, naturalFrequency, extremeLoading, pileDriving, suctionBucketInstallation, vesselImpact, reliabilityLifeExtension,
  ],
};
