// Suite 10 — Mooring Systems and Anchors.
import staticEquilibrium from './static-equilibrium.js';
import lineDynamics from './line-dynamics.js';
import coupled from './coupled-vessel-mooring.js';
import tensionOffset from './tension-offset.js';
import extremeEnvironment from './extreme-environment.js';
import lineFatigue from './line-fatigue.js';
import lineStrength from './line-strength.js';
import seabedContact from './seabed-contact.js';
import anchorHolding from './anchor-holding.js';
import dragAnchor from './drag-anchor.js';
import suctionAnchor from './suction-anchor.js';
import anchorPile from './anchor-pile.js';
import lineBreak from './line-break.js';
import snapLoads from './snap-loads.js';
import installation from './installation-pretension.js';
import clearance from './interference-clearance.js';
import reliability from './reliability.js';

export default {
  n: 10,
  modules: [staticEquilibrium, lineDynamics, coupled, tensionOffset, extremeEnvironment, lineFatigue, lineStrength, seabedContact, anchorHolding,
    dragAnchor, suctionAnchor, anchorPile, lineBreak, snapLoads, installation, clearance, reliability],
};
