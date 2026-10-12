// Suite 2 — Offshore topsides and deck structures.
import inplaceStrength from './inplace-strength.js';
import localFeStress from './local-fe-stress.js';
import deckPlateFraming from './deck-plate-framing.js';
import deflectionServiceability from './deflection-serviceability.js';
import bucklingStability from './buckling-stability.js';
import fatigueVibration from './fatigue-vibration.js';
import machineryFoundation from './machinery-foundation.js';
import cranePedestal from './crane-pedestal.js';
import skidGrillage from './skid-grillage.js';
import blastResponse from './blast-response.js';
import fireResponse from './fire-response.js';
import droppedObject from './dropped-object.js';
import windLoad from './wind-load.js';
import seismicResponse from './seismic-response.js';
import jacketInterface from './jacket-interface.js';
import liftInstallation from './lift-installation.js';
import seaTransport from './sea-transport.js';
import lifeExtension from './life-extension.js';
import reserveStrength from './reserve-strength.js';

export default {
  n: 2,
  modules: [
    inplaceStrength, localFeStress, deckPlateFraming, deflectionServiceability, bucklingStability, fatigueVibration, machineryFoundation, cranePedestal, skidGrillage,
    blastResponse, fireResponse, droppedObject, windLoad, seismicResponse, jacketInterface, liftInstallation, seaTransport, lifeExtension, reserveStrength,
  ],
};
