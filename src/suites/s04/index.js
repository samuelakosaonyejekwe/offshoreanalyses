// Suite 4 — Subsea support structures, mudmats and templates.
import onBottomStability from './on-bottom-stability.js';
import inplaceStrength from './inplace-strength.js';
import bearingCapacity from './bearing-capacity.js';
import slidingOverturning from './sliding-overturning.js';
import mudmatSoilInteraction from './mudmat-soil-interaction.js';
import settlement from './settlement.js';
import bucklingLocalStress from './buckling-local-stress.js';
import hydrodynamicLoading from './hydrodynamic-loading.js';
import skirtSuctionUplift from './skirt-suction-uplift.js';
import equipmentLoadDistribution from './equipment-load-distribution.js';
import pipelineJumperReactions from './pipeline-jumper-reactions.js';
import scour from './scour.js';
import loweringLanding from './lowering-landing.js';
import fatigueIntegrity from './fatigue-integrity.js';
import droppedObject from './dropped-object.js';

export default {
  n: 4,
  modules: [
    onBottomStability, inplaceStrength, bearingCapacity, slidingOverturning, mudmatSoilInteraction, settlement, bucklingLocalStress,
    hydrodynamicLoading, skirtSuctionUplift, equipmentLoadDistribution, pipelineJumperReactions, scour, loweringLanding, fatigueIntegrity, droppedObject,
  ],
};
