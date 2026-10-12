import waveTransformation from './wave-transformation.js';
import runupOvertopping from './runup-overtopping.js';
import armourStability from './armour-stability.js';
import wallStructure from './wall-structure.js';
import caissonWaveLoad from './caisson-wave-load.js';
import pierWaveLoad from './pier-wave-load.js';
import seawallStability from './seawall-stability.js';
import foundationSsi from './foundation-ssi.js';
import scourSediment from './scour-sediment.js';
import morphodynamics from './morphodynamics.js';
import seismicResponse from './seismic-response.js';
import fatigueDeterioration from './fatigue-deterioration.js';
import upliftFlotation from './uplift-flotation.js';
import floodInundation from './flood-inundation.js';
import intakeOutfall from './intake-outfall.js';
import vesselImpact from './vessel-impact.js';
import reliabilityServiceLife from './reliability-service-life.js';

export default {
  n: 7,
  modules: [waveTransformation, runupOvertopping, armourStability, wallStructure, caissonWaveLoad, pierWaveLoad, seawallStability, foundationSsi,
    scourSediment, morphodynamics, seismicResponse, fatigueDeterioration, upliftFlotation, floodInundation, intakeOutfall, vesselImpact, reliabilityServiceLife],
};
