// Suite 15 — Offshore installation and marine operations.
import craneLift from './crane-lift.js';
import rigging from './rigging.js';
import jacketLaunch from './jacket-launch.js';
import jacketUpending from './jacket-upending.js';
import loadout from './loadout.js';
import transportSeafastening from './transport-seafastening.js';
import bargeStability from './barge-stability.js';
import towing from './towing.js';
import pipeCableLay from './pipe-cable-lay.js';
import subseaLowering from './subsea-lowering.js';
import splashZone from './splash-zone.js';
import vesselMotion from './vessel-motion.js';
import pileDriving from './pile-driving.js';
import suctionInstallation from './suction-installation.js';
import floatover from './floatover.js';
import weatherWindow from './weather-window.js';
import droppedObjectCollision from './dropped-object-collision.js';
import temporaryStability from './temporary-stability.js';
import decommissioning from './decommissioning.js';

export default {
  n: 15,
  modules: [
    craneLift, rigging, jacketLaunch, jacketUpending, loadout, transportSeafastening, bargeStability, towing, pipeCableLay, subseaLowering,
    splashZone, vesselMotion, pileDriving, suctionInstallation, floatover, weatherWindow, droppedObjectCollision, temporaryStability, decommissioning,
  ],
};
