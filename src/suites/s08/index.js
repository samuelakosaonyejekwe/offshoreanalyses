// Suite 8 — Subsea Pipelines and Flowlines: mechanical design, pipe–soil and stability, installation, flow assurance
// and integrity modules.
import onBottomStability from './on-bottom-stability.js';
import hydrodynamicLoads from './hydrodynamic-loads.js';
import wallThickness from './wall-thickness.js';
import stress from './stress.js';
import localBuckling from './local-buckling.js';
import globalBuckling from './global-buckling.js';
import upheavalEmbedment from './upheaval-embedment.js';
import pipeSoil from './pipe-soil.js';
import thermalExpansion from './thermal-expansion.js';
import expansionWalking from './expansion-walking.js';
import freeSpan from './free-span.js';
import crossing from './crossing.js';
import installation from './installation.js';
import fatigueFracture from './fatigue-fracture.js';
import fluidPvt from './fluid-pvt.js';
import multiphaseFlow from './multiphase-flow.js';
import transientFlow from './transient-flow.js';
import slugging from './slugging.js';
import hydrate from './hydrate.js';
import waxScale from './wax-scale.js';
import corrosionErosion from './corrosion-erosion.js';
import integrity from './integrity.js';

export default {
  n: 8,
  modules: [wallThickness, stress, localBuckling, onBottomStability, hydrodynamicLoads, pipeSoil, globalBuckling, upheavalEmbedment, thermalExpansion, expansionWalking,
    freeSpan, crossing, installation, fatigueFracture, fluidPvt, multiphaseFlow, transientFlow, slugging, hydrate, waxScale, corrosionErosion, integrity],
};
