// Suite 5 — Wharves, Jetties, Quays and Marine Terminals.
import globalStrength from './global-strength.js';
import berthingEnergy from './berthing-energy.js';
import mooringBollard from './mooring-bollard.js';
import environmentalLoading from './environmental-loading.js';
import pileCapacity from './pile-capacity.js';
import quayWall from './quay-wall.js';
import soilStructure from './soil-structure.js';
import deckSlab from './deck-slab.js';
import seismicLiquefaction from './seismic-liquefaction.js';
import scour from './scour.js';
import deflectionSettlement from './deflection-settlement.js';
import fatigueVibration from './fatigue-vibration.js';
import limitStates from './limit-states.js';
import vesselOperational from './vessel-operational.js';
import accidentalCollapse from './accidental-collapse.js';
import deteriorationLife from './deterioration-life.js';

export default {
  n: 5,
  modules: [globalStrength, berthingEnergy, mooringBollard, environmentalLoading, pileCapacity, quayWall, soilStructure, deckSlab, seismicLiquefaction, scour,
    deflectionSettlement, fatigueVibration, limitStates, vesselOperational, accidentalCollapse, deteriorationLife],
};
