// Optical-fibre communication link in a submarine telecom cable or in the fibre element of a power cable: power
// budget with splice, connector, repair, ageing and bend losses, chromatic- and polarisation-mode-dispersion limits,
// single-mode waveguide parameters, nonlinear launch-power threshold, repeater count and power-feed voltage for
// repeatered systems, fibre strain life, and the reach of distributed temperature sensing on the same fibre.
import { linkBudget, spectralAttenuation, chromaticDispersion, dispersionLimitedLength, pmdLimitedLength, dispersionCompensation, stepIndexFibre, lp01Mode, macroBendLoss, osnr, repeaterSpacing, qFactorBer, nonlinearThresholds, powerFeedVoltage, sensingRange, fibreStaticFatigue } from '../../engine/electrical/optical.js';
import { indexProfile } from './_optical.js';
import { polylineLength, adv, check, worst, sev, fmt, round, run, need, linspace, CONFIRM, CO2_COPPER, CO2_MGO, MGO_T_PER_DAY } from './_common.js';

const E_GLASS = 72e9;         // Pa, Young's modulus of silica fibre (strain ↔ stress)
const YEAR = 31557600;

export default {
  id: 's13-optical-link', suite: 13, title: 'Optical link budget and fibre performance',
  summary: 'Attenuation, power budget, dispersion margins and maximum span of the fibre link, with repeater and power-feed sizing for long systems and the reach of fibre temperature monitoring.',
  covers: [9], fidelity: 'analytical',
  standards: ['ITU-T G.652 / G.654 / G.655 / G.657 (fibre attributes)', 'ITU-T G.957 / G.691 (optical interfaces and power budgets)', 'ITU-T G.973 (unrepeatered submarine systems)', 'ITU-T G.977 (optically amplified submarine systems, power feeding)', 'ITU-T G.972 (definitions for submarine cable systems)', 'IEC TR 62048 (fibre reliability, power-law theory)', 'Agrawal, Fiber-Optic Communication Systems (4th ed.) ch. 2, 4, 7'],
  models: ['Optical link budget (attenuation, splice, connector, repair and ageing losses)', 'Rayleigh scattering, infrared absorption and OH-peak attenuation spectrum', 'Chromatic dispersion (three-term Sellmeier fit) and dispersion-limited length', 'Polarisation-mode dispersion', 'Optical waveguide mode theory: step-index V number, cut-off, LP01 mode (Maxwell-based weakly guiding solution)', 'Refractive-index profile model', 'Macro-bend loss (Marcuse)', 'OSNR of an amplified chain and repeater spacing', 'Q factor and bit-error ratio', 'Stimulated Brillouin / Raman thresholds', 'Constant-current power feed of repeaters (Ohm\'s law on the cable conductor)', 'Power-law static fatigue of proof-tested fibre', 'Distributed sensing reach'],
  limits: 'Single channel, power-budget and bit-rate–distance rules: Kerr nonlinearity (self- and cross-phase modulation), WDM channel interaction, gain tilt, coherent-receiver digital compensation and forward-error-correction performance need a propagation simulation and vendor data and are outside the module. The bend-loss formula treats an infinite cladding and is for trends — use manufacturer macro-bend data for budgets. Fibre and equipment values are overridable typical figures to confirm against the supplier and the governing ITU-T edition.',
  publishes: ['cable.opticalMargin', 'cable.opticalMaxSpan'],
  searchTerms: 'submarine fibre optic link budget unrepeatered span dispersion OSNR power feed',
  inputs: [
    { key: 'route', label: 'Imported route (optional)', type: 'geometry', kinds: ['polyline'], default: null, group: 'Link', role: 'geometry', help: 'Route polyline; its length replaces the link length. Leave empty to use the form value.' },
    { key: 'L', label: 'Link length', unit: 'km', type: 'number', default: 60, min: 0.1, max: 15000, group: 'Link', role: 'geometry', bind: 'cable.length', design: true },
    { key: 'mode', label: 'System', type: 'select', default: 'unrepeatered', group: 'Link', role: 'bc', options: [{ value: 'unrepeatered', label: 'Unrepeatered (terminal equipment only)' }, { value: 'repeatered', label: 'Repeatered (optical amplifiers in the cable)' }] },
    { key: 'slack', label: 'Cable slack and fibre over-length', unit: '%', type: 'number', default: 2, min: 0, max: 15, group: 'Link', role: 'geometry', help: 'Fibre length beyond the route length: lay slack, helical lay of the fibre unit, jointing loops.' },
    { key: 'depth', label: 'Water depth for repairs', unit: 'm', type: 'number', default: 40, min: 1, max: 8000, group: 'Link', role: 'environment', bind: 'site.depth' },
    { key: 'lambda', label: 'Signal wavelength', unit: 'nm', type: 'number', default: 1550, min: 1260, max: 1650, group: 'Transmitter and receiver', role: 'ic', help: 'Initial optical operating condition of the link.' },
    { key: 'Ptx', label: 'Launched power per channel', unit: 'dBm', type: 'number', default: 3, min: -10, max: 30, group: 'Transmitter and receiver', role: 'ic', design: true, help: 'Transmitter output into the line fibre (after any booster).' },
    { key: 'sens', label: 'Receiver sensitivity at the required BER', unit: 'dBm', type: 'number', default: -24, min: -50, max: 0, group: 'Transmitter and receiver', role: 'bc', help: 'Receiver requirement: about −18 to −28 dBm at 10 Gbit/s depending on receiver type and FEC.' },
    { key: 'bitRate', label: 'Bit rate', unit: 'Gbit/s', type: 'number', default: 10, min: 0.1, max: 400, group: 'Transmitter and receiver', role: 'bc' },
    { key: 'srcWidth', label: 'Source spectral width (rms)', unit: 'nm', type: 'number', default: 0, min: 0, max: 5, group: 'Transmitter and receiver', role: 'ic', help: '0 = narrow, externally modulated laser (transform-limited limit applies).' },
    { key: 'gain', label: 'Remote / Raman amplification gain', unit: 'dB', type: 'number', default: 0, min: 0, max: 40, group: 'Transmitter and receiver', role: 'bc', showIf: { key: 'mode', in: ['unrepeatered'] } },
    { key: 'berReq', label: 'Required bit-error ratio', unit: '-', type: 'number', default: 1e-12, min: 1e-15, max: 1e-3, group: 'Transmitter and receiver', role: 'criteria' },
    { key: 'alpha', label: 'Cabled fibre attenuation', unit: 'dB/km', type: 'number', default: 0.2, min: 0.14, max: 1.0, group: 'Fibre', role: 'material', uncertain: { dist: 'normal', cov: 0.05 }, help: 'Typical cabled values at 1550 nm: G.652.D 0.20, G.654 0.17, G.655 0.22; 0.35 at 1310 nm.' },
    { key: 'rayleigh', label: 'Rayleigh scattering coefficient', unit: 'dB·µm⁴/km', type: 'number', default: 0.9, min: 0.6, max: 1.3, group: 'Fibre', role: 'material', help: 'For the spectral model: 0.75–1.0 for germanium-doped cores, lower for pure-silica cores.' },
    { key: 'lambda0', label: 'Zero-dispersion wavelength', unit: 'nm', type: 'number', default: 1310, min: 1260, max: 1600, group: 'Fibre', role: 'material' },
    { key: 'S0', label: 'Dispersion slope at the zero-dispersion wavelength', unit: 'ps/(nm²·km)', type: 'number', default: 0.092, min: 0.03, max: 0.12, group: 'Fibre', role: 'material' },
    { key: 'dcRatio', label: 'Dispersion compensated', unit: '%', type: 'number', default: 80, min: 0, max: 99.99, group: 'Fibre', role: 'bc', help: 'Share of the accumulated dispersion removed by compensating fibre, gratings or receiver processing.' },
    { key: 'pmd', label: 'PMD coefficient', unit: 'ps/√km', type: 'number', default: 0.1, min: 0.01, max: 2, group: 'Fibre', role: 'material', help: 'G.652.D link design value ≤ 0.20.' },
    { key: 'coreR', label: 'Core radius', unit: 'µm', type: 'number', default: 4.1, min: 2, max: 35, group: 'Fibre', role: 'geometry' },
    { key: 'n1', label: 'Core refractive index', unit: '-', type: 'number', default: 1.4504, min: 1.4, max: 1.6, group: 'Fibre', role: 'material' },
    { key: 'n2', label: 'Cladding refractive index', unit: '-', type: 'number', default: 1.4447, min: 1.39, max: 1.6, group: 'Fibre', role: 'material' },
    { key: 'Aeff', label: 'Effective area', unit: 'µm²', type: 'number', default: 80, min: 30, max: 160, group: 'Fibre', role: 'material', help: 'About 80 for G.652, 110–150 for G.654.' },
    { key: 'spliceLoss', label: 'Loss per splice', unit: 'dB', type: 'number', default: 0.05, min: 0, max: 0.5, group: 'Losses', role: 'bc' },
    { key: 'spliceInt', label: 'Section length between splices', unit: 'km', type: 'number', default: 20, min: 0.5, max: 200, group: 'Losses', role: 'geometry' },
    { key: 'connLoss', label: 'Loss per connector', unit: 'dB', type: 'number', default: 0.5, min: 0, max: 1.5, group: 'Losses', role: 'bc' },
    { key: 'nConn', label: 'Connectors', type: 'int', default: 2, min: 0, max: 20, group: 'Losses', role: 'bc' },
    { key: 'nRepairs', label: 'Repairs allowed for over the life', type: 'int', default: 3, min: 0, max: 50, group: 'Losses', role: 'operation', help: 'Each repair adds two splices and about 2.5 water depths of cable.' },
    { key: 'ageing', label: 'Ageing allowance', unit: 'dB/km', type: 'number', default: 0.005, min: 0, max: 0.05, group: 'Losses', role: 'material', help: 'Hydrogen and radiation-induced loss over the design life.' },
    { key: 'bendR', label: 'Smallest fibre bend radius in joints and terminations', unit: 'mm', type: 'number', default: 30, min: 5, max: 200, group: 'Losses', role: 'geometry' },
    { key: 'turns', label: 'Fibre turns stored at that radius', type: 'int', default: 20, min: 0, max: 2000, group: 'Losses', role: 'geometry' },
    { key: 'sysMargin', label: 'System margin', unit: 'dB', type: 'number', default: 3, min: 0, max: 10, group: 'Acceptance criteria', role: 'criteria', help: `End-of-life operating margin kept in hand (ITU-T G.957 practice 3 dB). ${CONFIRM}` },
    { key: 'dgdFrac', label: 'Allowed DGD as a fraction of the bit period', unit: '-', type: 'number', default: 0.1, min: 0.02, max: 0.5, group: 'Acceptance criteria', role: 'criteria' },
    { key: 'sbsMargin', label: 'Brillouin threshold increase from modulation / dithering', unit: 'dB', type: 'number', default: 10, min: 0, max: 30, group: 'Acceptance criteria', role: 'criteria', help: 'The threshold formula is for an unmodulated narrow line; modulated signals and dithered sources tolerate more.' },
    { key: 'nf', label: 'Amplifier noise figure', unit: 'dB', type: 'number', default: 4.5, min: 3, max: 10, group: 'Repeatered system', role: 'material', showIf: { key: 'mode', in: ['repeatered'] } },
    { key: 'osnrReq', label: 'Required OSNR at the receiver', unit: 'dB', type: 'number', default: 15, min: 5, max: 35, group: 'Repeatered system', role: 'criteria', showIf: { key: 'mode', in: ['repeatered'] } },
    { key: 'vRep', label: 'Voltage drop per repeater', unit: 'V', type: 'number', default: 40, min: 5, max: 200, group: 'Repeatered system', role: 'bc', showIf: { key: 'mode', in: ['repeatered'] } },
    { key: 'iFeed', label: 'Line current', unit: 'A', type: 'number', default: 1.0, min: 0.3, max: 2, group: 'Repeatered system', role: 'bc', showIf: { key: 'mode', in: ['repeatered'] } },
    { key: 'rCond', label: 'Power-conductor resistance', unit: 'Ω/km', type: 'number', default: 1.0, min: 0.3, max: 3, group: 'Repeatered system', role: 'material', showIf: { key: 'mode', in: ['repeatered'] } },
    { key: 'pfeMax', label: 'Power-feed equipment voltage rating per station', unit: 'kV', type: 'number', default: 15, min: 1, max: 20, group: 'Repeatered system', role: 'criteria', showIf: { key: 'mode', in: ['repeatered'] } },
    { key: 'tension', label: 'Sustained cable tension', unit: 'kN', type: 'number', default: 30, min: 0, max: 2000, group: 'Fibre strain', role: 'bc', help: 'Permanent (residual) tension after laying; suspended sections carry their own weight.' },
    { key: 'EA', label: 'Cable axial stiffness', unit: 'MN', type: 'number', default: 40, min: 1, max: 5000, group: 'Fibre strain', role: 'material', help: 'Lightweight telecom cable 20–50 MN; armoured power cable several hundred.' },
    { key: 'efl', label: 'Fibre strain-free window (excess fibre length)', unit: '%', type: 'number', default: 0.05, min: 0, max: 1, group: 'Fibre strain', role: 'geometry', help: 'Cable elongation taken up before the fibres strain: near zero in tight structures, 0.2–0.5 % in loose tubes.' },
    { key: 'proof', label: 'Fibre proof-test strain', unit: '%', type: 'number', default: 1.0, min: 0.5, max: 2.5, group: 'Fibre strain', role: 'material', help: '1 % (0.69 GPa) is standard; 2 % for submarine-grade fibre.' },
    { key: 'nFat', label: 'Stress-corrosion parameter n', unit: '-', type: 'number', default: 20, min: 15, max: 30, group: 'Fibre strain', role: 'material' },
    { key: 'life', label: 'Design life', unit: 'years', type: 'number', default: 25, min: 1, max: 50, group: 'Fibre strain', role: 'criteria', bind: 'life.designYears' },
    { key: 'dts', label: 'Distributed temperature sensing on a spare fibre', type: 'bool', default: true, group: 'Monitoring', role: 'operation' },
    { key: 'dtsRange', label: 'Interrogator dynamic range (two-way)', unit: 'dB', type: 'number', default: 30, min: 5, max: 60, group: 'Monitoring', role: 'material', showIf: { key: 'dts', in: [true] } },
    { key: 'dtsLoop', label: 'Interrogated from both ends', type: 'bool', default: false, group: 'Monitoring', role: 'operation', showIf: { key: 'dts', in: [true] } },
  ],
  compute(inp) {
    const Lroute = polylineLength(inp.route) ?? inp.L * 1e3, Lf = Lroute * (1 + inp.slack / 100);   // fibre length, m
    need(inp.n1 > inp.n2, 'The core refractive index must exceed the cladding index: no guided mode otherwise.');
    const B = inp.bitRate * 1e9, rep = inp.mode === 'repeatered';
    // waveguide
    const wg = stepIndexFibre({ coreRadius: inp.coreR * 1e-6, n1: inp.n1, n2: inp.n2, wavelengthNm: inp.lambda });
    const lp = lp01Mode(wg.V);
    const bend = run(() => macroBendLoss({ bendRadius: inp.bendR / 1e3, coreRadius: inp.coreR * 1e-6, n1: inp.n1, n2: inp.n2, wavelengthNm: inp.lambda, turns: inp.turns }), 'Macro-bend loss');
    const bendDb = Number.isFinite(bend.lossDb) ? bend.lossDb : 0;
    // losses and budget
    const repairDb = inp.nRepairs * (2 * inp.spliceLoss + 2.5 * inp.depth / 1e3 * inp.alpha), ageDb = inp.ageing * Lf / 1e3;
    const budgetArgs = { attenuationDbPerKm: inp.alpha, launchPowerDbm: inp.Ptx, sensitivityDbm: inp.sens, spliceLossDb: inp.spliceLoss, spliceInterval: inp.spliceInt * 1e3, connectorLossDb: inp.connLoss, nConnectors: inp.nConn, repairMarginDb: repairDb, systemMarginDb: inp.sysMargin, otherLossDb: bendDb, amplifierGainDb: rep ? 0 : inp.gain };
    const lb = run(() => linkBudget({ ...budgetArgs, length: Lf, ageingDb: ageDb }), 'Link budget');
    // the ageing allowance grows with length: maximum span by bisection on the margin
    const marginAt = (len) => linkBudget({ ...budgetArgs, length: len, ageingDb: inp.ageing * len / 1e3 }).marginDb;
    let maxSpan = 0;
    if (marginAt(1) > 0) { let lo = 1, hi = 2e3; while (marginAt(hi) > 0 && hi < 5e7) hi *= 2; for (let i = 0; i < 80; i++) { const mid = 0.5 * (lo + hi); if (marginAt(mid) >= 0) lo = mid; else hi = mid; } maxSpan = lo / (1 + inp.slack / 100); }
    // spectrum
    const spec = spectralAttenuation({ wavelengthNm: inp.lambda, rayleigh: inp.rayleigh });
    // dispersion
    const cd = chromaticDispersion({ wavelengthNm: inp.lambda, zeroDispersionNm: inp.lambda0, slopePsPerNm2Km: inp.S0 });
    const resid = 1 - inp.dcRatio / 100;
    const dl = dispersionLimitedLength({ bitRate: B, DPsPerNmKm: cd.DPsPerNmKm, wavelengthNm: inp.lambda, ...(inp.srcWidth > 0 ? { sourceWidthNm: inp.srcWidth } : {}) });
    const Ldisp = Number.isFinite(dl.length) ? dl.length / resid : 1e12;
    const accum = cd.DPsPerNmKm * Lf / 1e3;
    const dcf = inp.dcRatio > 0 && Math.abs(cd.DPsPerNmKm) > 1e-6 ? dispersionCompensation({ length: Lf, DPsPerNmKm: cd.DPsPerNmKm, DdcfPsPerNmKm: cd.DPsPerNmKm > 0 ? -100 : 100, ratio: inp.dcRatio / 100 }) : { lengthDcf: 0, residualPsPerNm: accum };
    const pm = pmdLimitedLength({ bitRate: B, pmdPsPerSqrtKm: inp.pmd, fraction: inp.dgdFrac, length: Lf });
    // nonlinear threshold (one span for a repeatered system)
    const q = qFactorBer({ ber: inp.berReq });
    let rp = null, pfe = null, osnrDb = null;
    if (rep) {
      rp = run(() => repeaterSpacing({ length: Lf, attenuationDbPerKm: inp.alpha, channelPowerDbm: inp.Ptx, noiseFigureDb: inp.nf, requiredOsnrDb: inp.osnrReq, extraSpanLossDb: 2 * inp.spliceLoss + inp.ageing * 80 }), 'No repeater spacing meets the required OSNR — raise the launch power or lower the requirement');
      osnrDb = osnr({ channelPowerDbm: inp.Ptx, noiseFigureDb: inp.nf, spanLossDb: rp.spanLossDb, nSpans: rp.nSpans });
      pfe = powerFeedVoltage({ nRepeaters: rp.nRepeaters, repeaterVoltage: inp.vRep, current: inp.iFeed, conductorResistance: inp.rCond / 1e3, length: Lroute, earthPotentialMargin: 0.1 * Lroute / 1e3, feed: 'double' });
    }
    const nl = nonlinearThresholds({ length: rep ? rp.spacing : Lf, attenuationDbPerKm: inp.alpha, effectiveArea: inp.Aeff * 1e-12 });
    const sbsLim = nl.sbsThresholdDbm + inp.sbsMargin;
    // fibre strain and static fatigue
    const cableStrain = inp.tension * 1e3 / (inp.EA * 1e6), fibreStrain = Math.max(0, cableStrain - inp.efl / 100);
    const fat = fibreStaticFatigue({ proofStress: E_GLASS * inp.proof / 100, proofTime: 1, n: inp.nFat, designLife: inp.life * YEAR, ...(fibreStrain > 0 ? { serviceStress: E_GLASS * fibreStrain } : {}) });
    const strainAllow = fat.allowableStress / E_GLASS;
    // monitoring reach
    const reach = inp.dts ? sensingRange({ dynamicRangeDb: inp.dtsRange, attenuationDbPerKm: inp.alpha, passes: 2, fixedLossDb: inp.nConn * inp.connLoss / 2 }) * (inp.dtsLoop ? 2 : 1) : 0;
    const checks = [
      ...(!rep ? [check('Optical power budget', lb.totalLossDb + inp.sysMargin, Math.max(inp.Ptx + (inp.gain || 0) - inp.sens, 1e-9), 'dB', 'ITU-T G.957 Appendix / G.973 (budget method)', 'End-of-life link loss plus the system margin against the gain available between transmitter and receiver. When it is exceeded the received power falls below the receiver sensitivity and the bit-error ratio rises.')] : [
        check('Optical signal-to-noise ratio', inp.osnrReq, osnrDb, 'dB', 'ITU-T G.977; OSNR = 58 + P − NF − Lspan − 10·log N', 'Required OSNR against the OSNR delivered by the amplifier chain. Noise from every amplifier accumulates; more, shorter spans are the remedy.'),
        check('Power-feed voltage per station', pfe.voltagePerStation / 1e3, inp.pfeMax, 'kV', 'ITU-T G.977 (power feeding); Ohm\'s law on the cable conductor', 'Voltage each terminal station must apply to drive the line current through repeaters and conductor against the feed-equipment rating and the cable insulation.'),
      ]),
      check('Chromatic dispersion limit', Lf / 1e3, Ldisp / 1e3, 'km', `Agrawal §2.4.3 (${dl.criterion}); ITU-T G.691`, 'Fibre length against the length at which pulse spreading closes the eye, after the stated compensation. Beyond it the receiver cannot separate neighbouring bits whatever the power.'),
      check('Polarisation-mode dispersion', pm.dgdPs, inp.dgdFrac * 1e12 / B, 'ps', 'ITU-T G.691 Appendix I; Agrawal §2.4.5', 'Mean differential group delay of the link against the allowed fraction of the bit period. PMD is random and cannot be compensated passively.'),
      check('Launch power against the Brillouin threshold', 10 ** (inp.Ptx / 10), 10 ** (sbsLim / 10), 'mW', 'Smith (1972); Agrawal §2.6', 'Launched power against the stimulated-Brillouin threshold with the modulation allowance. Above it power is scattered back and the received power saturates.'),
      check('Single-mode operation', wg.cutoffWavelengthNm, inp.lambda, 'nm', 'ITU-T G.650.1 (cut-off wavelength); V < 2.405', 'Theoretical cut-off wavelength against the signal wavelength. Below cut-off a second mode propagates and modal dispersion destroys the bandwidth.'),
      check('Sustained fibre strain', fibreStrain * 100 + 1e-9, strainAllow * 100, '%', 'IEC TR 62048 (power-law static fatigue); ITU-T G.972 (NPTS)', 'Permanent tensile strain of the fibres against the strain a proof-tested fibre survives for the design life. Glass under sustained stress in a humid environment fails by slow crack growth.'),
      ...(inp.dts ? [check('Temperature-sensing reach', Lf / 1e3, Math.max(reach / 1e3, 1e-9), 'km', 'Interrogator optical budget (two-way loss); IEC 61757-2-2 (DTS)', 'Fibre length against the distance the interrogator can measure. Beyond the reach the far part of the cable — often the landfall hotspot — is not monitored.')] : []),
    ];
    // the waveguide check is a property of the fibre, not a utilisation that grows with the link: it governs only when it fails
    const live = checks.filter((c) => c.uc > 1 || !c.name.startsWith('Single-mode'));
    const gov = worst(live.length ? live : checks);
    // power level along the link (fibre + evenly spread splices), end-of-life
    const xs = linspace(0, Lf / 1e3, 61), perKm = (lb.totalLossDb - inp.nConn * inp.connLoss) / (Lf / 1e3);
    const level = xs.map((x) => inp.Ptx + (rep ? 0 : inp.gain) - inp.nConn * inp.connLoss / 2 - perKm * x - (x >= Lf / 1e3 ? inp.nConn * inp.connLoss / 2 : 0));
    const lams = linspace(1250, 1650, 41), radii = linspace(Math.max(5, inp.bendR * 0.3), Math.max(60, inp.bendR * 1.5), 30);
    const rr = linspace(0, 3 * inp.coreR, 40), dn = wg.delta;
    const Ls = linspace(Math.max(1, 0.2 * Lf / 1e3), 2 * Lf / 1e3, 30);
    const vesselDaysRepeater = rep ? rp.nRepeaters * 0.1 : 0;
    return {
      kpis: [
        { key: 'fibreLoss', label: 'Fibre attenuation of the link', value: round(lb.fibreLossDb, 5), unit: 'dB', explain: 'Attenuation coefficient × fibre length including slack: the dominant loss of a long link (end to end; in a repeatered system the amplifiers restore it span by span).' },
        { key: 'totalLoss', label: 'End-of-life link loss', value: round(lb.totalLossDb, 5), unit: 'dB', explain: 'Fibre, splices, connectors, repair allowance, ageing and bend losses added: what the transmitter and receiver must bridge.' },
        { key: 'Prx', label: 'Received power', value: round(lb.receivedPowerDbm, 5), unit: 'dBm', explain: 'Launched power less the link loss (plus any remote gain); compare with the receiver sensitivity.' },
        { key: 'margin', label: 'Power margin after the system margin', value: round(lb.marginDb, 5), unit: 'dB', ...(rep ? {} : { limit: 0, sense: 'min' }), explain: 'Received power minus sensitivity minus the system margin. Positive means the budget closes with the margin intact; for a repeatered system it is shown for one unamplified pass for reference only.' },
        { key: 'maxSpan', label: 'Maximum unrepeatered route length', value: round(maxSpan / 1e3), unit: 'km', explain: 'Route length at which the power margin falls to zero with these losses: beyond it the system needs remote amplification or repeaters.' },
        { key: 'nSplices', label: 'Splices in the link', value: lb.nSplices, unit: '-', explain: 'Section splices plus two terminal splices; repairs add two each within the repair allowance.' },
        { key: 'alphaModel', label: 'Attenuation from the spectral model', value: round(spec.total, 5), unit: 'dB/km', explain: 'Rayleigh scattering plus infrared absorption at the signal wavelength: the intrinsic floor of silica fibre, for comparison with the cabled value entered.' },
        { key: 'rayleighLoss', label: 'Rayleigh scattering loss', value: round(spec.rayleigh, 5), unit: 'dB/km', explain: 'A/λ⁴ part of the attenuation; it falls steeply with wavelength, which is why long links use 1550 nm.' },
        { key: 'D', label: 'Chromatic dispersion', value: round(cd.DPsPerNmKm, 5), unit: 'ps/(nm·km)', explain: 'Spread of arrival time per nanometre of source width per kilometre at the signal wavelength.' },
        { key: 'dispAccum', label: 'Residual accumulated dispersion', value: round(dcf.residualPsPerNm), unit: 'ps/nm', explain: 'Dispersion left at the receiver after the stated compensation.' },
        { key: 'Ldisp', label: 'Dispersion-limited length', value: round(Math.min(Ldisp, 1e9) / 1e3), unit: 'km', explain: 'Fibre length at which chromatic dispersion limits the bit rate, with the compensation credited.' },
        { key: 'dgd', label: 'Mean differential group delay', value: round(pm.dgdPs, 5), unit: 'ps', explain: 'PMD coefficient × √length; it must stay a small fraction of the bit period.' },
        { key: 'V', label: 'Normalised frequency V', value: round(wg.V, 5), unit: '-', explain: 'Waveguide parameter 2π·a·NA/λ; the fibre is single-mode for V below 2.405.' },
        { key: 'cutoff', label: 'Theoretical cut-off wavelength', value: round(wg.cutoffWavelengthNm, 5), unit: 'nm', explain: 'Wavelength above which only the fundamental mode is guided.' },
        { key: 'mfd', label: 'Mode-field diameter', value: round(Number.isFinite(wg.modeFieldDiameter) ? wg.modeFieldDiameter * 1e6 : 0), unit: 'µm', explain: 'Gaussian width of the fundamental mode (Marcuse fit); mismatch between fibres causes splice loss. 0 when V is outside the fit range.' },
        { key: 'bLP01', label: 'Normalised propagation constant of LP01', value: round(lp.b), unit: '-', explain: 'Fraction of the index step seen by the fundamental mode; low values mean a weakly confined mode that is sensitive to bending.' },
        { key: 'bendLoss', label: 'Bend loss of the stored turns', value: round(bendDb, 4), unit: 'dB', explain: 'Macro-bend loss of the fibre loops in joints and terminations at the smallest radius (trend formula).' },
        { key: 'Qreq', label: 'Q factor for the required BER', value: round(q.Q, 5), unit: '-', explain: 'Signal-to-noise quality the receiver decision needs for the bit-error ratio: Q = 6 gives 10⁻⁹, Q = 7 gives 10⁻¹².' },
        { key: 'sbs', label: 'Brillouin threshold (unmodulated)', value: round(nl.sbsThresholdDbm), unit: 'dBm', explain: 'Launch power at which stimulated Brillouin back-scatter sets in for a narrow unmodulated line.' },
        { key: 'fibreStrain', label: 'Sustained fibre strain', value: round(fibreStrain * 100, 4), unit: '%', limit: round(strainAllow * 100, 4), sense: 'max', explain: 'Cable elongation at the sustained tension less the strain-free window of the fibre.' },
        { key: 'strainAllow', label: 'Allowable sustained fibre strain', value: round(strainAllow * 100, 5), unit: '%', explain: 'Strain a fibre proof-tested at the entered level survives for the design life (minimum-life form); commonly a fifth to a third of the proof strain.' },
        ...(inp.dts ? [{ key: 'dtsReach', label: 'Temperature-sensing reach', value: round(reach / 1e3), unit: 'km', explain: 'Fibre length over which the interrogator returns usable temperature data; double-ended interrogation covers twice the distance.' }] : []),
        ...(rep ? [
          { key: 'nRepeaters', label: 'Repeaters', value: rp.nRepeaters, unit: '-', explain: 'Smallest number of equally spaced optical amplifiers that meets the required OSNR.' },
          { key: 'spacing', label: 'Repeater spacing', value: round(rp.spacing / 1e3), unit: 'km', explain: 'Fibre length of one amplified span.' },
          { key: 'osnr', label: 'OSNR at the receiver', value: round(osnrDb), unit: 'dB', limit: inp.osnrReq, sense: 'min', explain: 'Optical signal-to-noise ratio in 0.1 nm after the amplifier chain.' },
          { key: 'pfeVoltage', label: 'Power-feed voltage per station', value: round(pfe.voltagePerStation / 1e3), unit: 'kV', limit: inp.pfeMax, sense: 'max', explain: 'Half of the system voltage with double-end feeding; it stresses the cable insulation continuously.' },
          { key: 'pfePower', label: 'Power fed into the cable', value: round(pfe.power / 1e3), unit: 'kW', explain: 'System voltage × line current: electrical loading of the power conductor, mostly dissipated in the conductor resistance.' },
        ] : []),
      ],
      checks,
      plots: [
        { type: 'line', title: 'Optical power along the link (end of life)', x: { label: 'Distance from the transmitter', unit: 'km' }, y: { label: 'Power', unit: 'dBm' },
          series: [{ name: 'Signal level', x: xs, y: level }, { name: 'Receiver sensitivity', x: [0, xs[60]], y: [inp.sens, inp.sens], dash: true }, { name: 'Sensitivity + system margin', x: [0, xs[60]], y: [inp.sens + inp.sysMargin, inp.sens + inp.sysMargin], dash: true }],
          explain: 'The level falls linearly in decibels with distance. For an unrepeatered link it must reach the far end above the upper dashed line; where it crosses that line is the maximum span.' },
        { type: 'line', title: 'Power margin against link length', x: { label: 'Fibre length', unit: 'km' }, y: { label: 'Margin', unit: 'dB' },
          series: [{ name: 'Margin', x: Ls, y: Ls.map((l) => marginAt(l * 1e3)) }, { name: 'Zero margin', x: [Ls[0], Ls[29]], y: [0, 0], dash: true }, { name: 'This link', x: [Lf / 1e3, Lf / 1e3], y: [Math.min(0, lb.marginDb), Math.max(0, lb.marginDb)], dash: true }],
          explain: 'Margin left as the link is lengthened. The crossing of zero is the reach of the terminal equipment; each 0.01 dB/km of attenuation moves it by several kilometres.' },
        { type: 'line', title: 'Attenuation spectrum of the fibre (model)', x: { label: 'Wavelength', unit: 'nm' }, y: { label: 'Attenuation', unit: 'dB/km' },
          series: [{ name: 'Total', x: lams, y: lams.map((l) => spectralAttenuation({ wavelengthNm: l, rayleigh: inp.rayleigh }).total) }, { name: 'Rayleigh scattering', x: lams, y: lams.map((l) => spectralAttenuation({ wavelengthNm: l, rayleigh: inp.rayleigh }).rayleigh), dash: true }, { name: 'Cabled value entered', x: [inp.lambda], y: [inp.alpha] }],
          explain: 'Rayleigh scattering falls with λ⁴ and infrared absorption rises beyond 1600 nm, leaving the low-loss window near 1550 nm. The cabled value sits above the model by the cabling and splice-free excess loss.' },
        { type: 'line', title: 'Bend loss against bend radius', x: { label: 'Bend radius', unit: 'mm' }, y: { label: 'log₁₀ loss per turn', unit: 'dB' },
          series: [{ name: 'Loss per turn', x: radii, y: radii.map((r) => Math.log10(Math.max(macroBendLoss({ bendRadius: r / 1e3, coreRadius: inp.coreR * 1e-6, n1: inp.n1, n2: inp.n2, wavelengthNm: inp.lambda }).lossDbPerTurn, 1e-30))) }, { name: 'Stored radius', x: [inp.bendR, inp.bendR], y: [-30, 1], dash: true }],
          explain: 'Bend loss rises exponentially as the radius shrinks: the knee marks the smallest safe storage radius for joint trays and for the cable minimum bend radius.' },
        { type: 'line', title: 'Refractive-index profile', x: { label: 'Radius', unit: 'µm' }, y: { label: 'Refractive index', unit: '-' },
          series: [{ name: 'Step profile', x: rr, y: rr.map((r) => indexProfile({ r: r * 1e-6, coreRadius: inp.coreR * 1e-6, n1: inp.n1, delta: Math.min(Math.max(dn, 1e-6), 0.49) })) }],
          explain: 'Index step between core and cladding that guides the light. The step height and core radius together set V, the cut-off wavelength and the bend sensitivity.' },
      ],
      tables: [
        { title: 'Loss budget', explain: 'One row per loss item of the end-of-life budget, in the order light meets them, with the resulting received power and margin.',
          columns: [{ key: 'item', label: 'Item' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ item: 'Launched power', v: inp.Ptx, u: 'dBm' }, ...(!rep && inp.gain > 0 ? [{ item: 'Remote / Raman gain', v: inp.gain, u: 'dB' }] : []), { item: `Fibre ${fmt(Lf / 1e3, 5)} km × ${inp.alpha} dB/km`, v: round(lb.fibreLossDb), u: 'dB' }, { item: `Splices (${lb.nSplices})`, v: round(lb.spliceLossDb), u: 'dB' }, { item: `Connectors (${inp.nConn})`, v: round(lb.connectorLossDb), u: 'dB' },
            { item: `Repair allowance (${inp.nRepairs} repairs)`, v: round(repairDb), u: 'dB' }, { item: 'Ageing', v: round(ageDb), u: 'dB' }, { item: 'Bends', v: round(bendDb, 3), u: 'dB' }, { item: 'Total loss', v: round(lb.totalLossDb), u: 'dB' }, { item: 'Received power', v: round(lb.receivedPowerDbm), u: 'dBm' },
            { item: 'Receiver sensitivity', v: inp.sens, u: 'dBm' }, { item: 'System margin', v: inp.sysMargin, u: 'dB' }, { item: 'Remaining margin', v: round(lb.marginDb), u: 'dB' }] },
        { title: 'Dispersion and waveguide summary', explain: 'Dispersion limits and single-mode parameters of the fibre at the signal wavelength, for comparison with the fibre data sheet.',
          columns: [{ key: 'q', label: 'Quantity' }, { key: 'v', label: 'Value' }, { key: 'u', label: 'Unit' }],
          rows: [{ q: 'Dispersion D', v: round(cd.DPsPerNmKm), u: 'ps/(nm·km)' }, { q: 'Accumulated dispersion', v: round(accum), u: 'ps/nm' }, { q: 'Compensating fibre length (−100 ps/(nm·km))', v: round(dcf.lengthDcf / 1e3), u: 'km' }, { q: 'Dispersion-limited length', v: round(Math.min(Ldisp, 1e9) / 1e3), u: 'km' }, { q: 'PMD-limited length', v: round(pm.length / 1e3), u: 'km' },
            { q: 'Numerical aperture', v: round(wg.NA), u: '-' }, { q: 'V number', v: round(wg.V), u: '-' }, { q: 'Cut-off wavelength', v: round(wg.cutoffWavelengthNm), u: 'nm' }, { q: 'Effective length (nonlinear)', v: round(nl.effectiveLength / 1e3), u: 'km' }, { q: 'Raman threshold', v: round(nl.srsThresholdDbm), u: 'dBm' }] },
      ],
      notes: [
        `Governing: ${gov.name} (utilisation ${fmt(gov.uc)}).`,
        rep ? `Amplified system: span loss ${fmt(rp.spanLossDb)} dB over ${fmt(rp.spacing / 1e3)} km restored by each repeater; the same terminal equipment would reach ${fmt(maxSpan / 1e3)} km unrepeatered.`
          : `Budget: ${fmt(lb.totalLossDb)} dB end-of-life loss over ${fmt(Lf / 1e3)} km of fibre; received ${fmt(lb.receivedPowerDbm)} dBm against ${inp.sens} dBm; margin ${fmt(lb.marginDb)} dB after the ${inp.sysMargin} dB system margin. Maximum unrepeatered route ${fmt(maxSpan / 1e3)} km.`,
        `Dispersion ${fmt(cd.DPsPerNmKm)} ps/(nm·km), limit ${fmt(Math.min(Ldisp, 1e9) / 1e3)} km with ${inp.dcRatio} % compensated (${dl.criterion}); DGD ${fmt(pm.dgdPs)} ps of ${fmt(inp.dgdFrac * 1e12 / B)} ps allowed.`,
        rep ? `Repeatered: ${rp.nRepeaters} repeaters at ${fmt(rp.spacing / 1e3)} km, OSNR ${fmt(osnrDb)} dB, ${fmt(pfe.voltagePerStation / 1e3)} kV per station at ${inp.iFeed} A.` : `Fibre V = ${fmt(wg.V)}, cut-off ${fmt(wg.cutoffWavelengthNm)} nm: ${wg.singleMode ? 'single-mode' : 'multi-mode'} at ${inp.lambda} nm.`,
      ],
      advice: [
        adv('engineering', sev(gov.uc), `${gov.name} governs at ${fmt(gov.uc)}. ${rep ? `${rp.nRepeaters} repeater(s) at ${fmt(rp.spacing / 1e3)} km give an OSNR of ${fmt(osnrDb)} dB (required ${inp.osnrReq} dB) with ${fmt(pfe.voltagePerStation / 1e3)} kV per station` : `Link loss ${fmt(lb.totalLossDb)} dB, margin ${fmt(lb.marginDb)} dB`}, dispersion limit ${fmt(Math.min(Ldisp, 1e9) / 1e3)} km for ${fmt(Lf / 1e3)} km, sustained fibre strain ${fmt(fibreStrain * 100)} % of ${fmt(strainAllow * 100)} % allowable.`,
          gov.uc > 1 ? (gov.name.startsWith('Optical power') ? `Close ${fmt(-lb.marginDb)} dB: lower-loss fibre (G.654 at 0.17 dB/km saves ${fmt((inp.alpha - 0.17) * Lf / 1e3)} dB here), more launch power within the nonlinear limit, a more sensitive receiver or remote amplification; otherwise repeaters are needed beyond ${fmt(maxSpan / 1e3)} km.` : gov.name.startsWith('Chromatic') ? `Compensate at least ${fmt((1 - dl.length / Lf) * 100, 5)} % of the dispersion, use a lower bit rate per channel or coherent detection with electronic compensation.` : gov.name.startsWith('Polarisation') ? 'Specify a lower PMD link design value for the fibre or reduce the bit rate per channel.' : gov.name.startsWith('Launch') ? `Reduce the launch power to ${fmt(sbsLim)} dBm or broaden the source linewidth (dithering).` : gov.name.startsWith('Single-mode') ? 'Operate above the cut-off wavelength or use a fibre with a smaller core / index step.' : gov.name.startsWith('Sustained') ? 'Reduce the residual lay tension, use a loose-tube design with more excess fibre length, or specify fibre proof-tested at 2 %.' : gov.name.startsWith('Temperature') ? 'Interrogate from both ends, use a longer-range interrogator or accept monitoring of the near part only.' : 'Shorten the amplifier spans or raise the feed-equipment rating; check the cable insulation for the higher voltage.')
            : 'Confirm the budget with the supplier\'s cabled attenuation and the as-built splice record, keep the repair allowance intact through commissioning, and re-check the margin at the end-of-life temperature and with the final equipment data.',
          'ITU-T G.957 / G.973 / G.977; ITU-T G.652–G.657; IEC TR 62048'),
        adv('economic', 'info', rep ? `The link needs ${rp.nRepeaters} repeater(s) and ${fmt(pfe.power / 1e3)} kW of line power (${fmt(pfe.power * 8.766 / 1e3)} MWh a year).` : `The link closes without repeaters with ${fmt(lb.marginDb)} dB in hand; the maximum unrepeatered route is ${fmt(maxSpan / 1e3)} km.`,
          rep ? `Each repeater is a submerged plant item to buy, qualify and — if it fails — recover: a 0.01 dB/km better fibre lengthens the spans by about ${fmt(0.01 * rp.spacing / 1e3 / inp.alpha)} km each and may remove ${Math.max(0, Math.round(rp.nRepeaters * 0.01 / inp.alpha))} repeater(s). Compare that with the fibre premium over the whole length.`
            : `Margin is an asset: ${fmt(lb.marginDb)} dB is worth ${fmt(lb.marginDb / (2 * inp.spliceLoss + 2.5 * inp.depth / 1e3 * inp.alpha))} further repairs before the link fails its budget, or an upgrade in bit rate later. Avoiding repeaters removes the power conductor, the feed equipment and their maintenance from the system.`,
          'Whole-life costing ISO 15686-5; value-of-information for spare margin and fibre count',
          rep ? `Line power ${fmt(pfe.power * 8.766 / 1e3)} MWh/yr; the power conductor adds about ${fmt(1.7241e-8 / (inp.rCond / 1e3) * 8890 * Lroute / 1e3)} t of copper ≈ ${fmt(1.7241e-8 / (inp.rCond / 1e3) * 8890 * Lroute / 1e3 * CO2_COPPER)} t CO₂e (screening ${CO2_COPPER} t CO₂e/t). An unrepeatered design avoids both.`
            : `No submerged electronics and no power conductor: the lowest material and energy footprint for the link. Each avoided fibre repair saves a cable-ship call of several days at about ${fmt(MGO_T_PER_DAY * CO2_MGO)} t CO₂ per day (screening: ${MGO_T_PER_DAY} t/day fuel, ${CO2_MGO} t CO₂ per t fuel).`),
      ],
      publish: { 'cable.opticalMargin': lb.marginDb, 'cable.opticalMaxSpan': maxSpan },
      quantities: { cable_km: round(Lroute / 1e3, 6), ...(rep ? { copper_t: round(1.7241e-8 / (inp.rCond / 1e3) * 8890 * Lroute / 1e3), vessel_days: { repeater_integration: round(vesselDaysRepeater) } } : {}) },
    };
  },
  verify: [
    // Budget by hand: 100 km of fibre (no slack) at 0.2 dB/km = 20 dB; sections of 5 km → 19 section splices + 2 terminal = 21 × 0.05 = 1.05 dB;
    // 2 connectors × 0.5 = 1 dB; no repairs, ageing or bends → total 22.05 dB. Received = 3 − 22.05 = −19.05 dBm;
    // margin = −19.05 − (−30) − 3 = 7.95 dB.
    { name: 'Optical power budget: loss, received power and margin',
      source: 'ITU-T G.957 Appendix (power budget method); hand calculation',
      inputs: { L: 100, slack: 0, spliceInt: 5, nRepairs: 0, ageing: 0, turns: 0, sens: -30, dcRatio: 90 },
      expect: [{ kpi: 'fibreLoss', value: 20, tol: 1e-6 }, { kpi: 'totalLoss', value: 22.05, tol: 1e-4 }, { kpi: 'Prx', value: -19.05, tol: 1e-4 }, { kpi: 'margin', value: 7.95, tol: 1e-3 }, { kpi: 'nSplices', value: 21, abs: 0 }] },
    // Fibre physics by hand at 1550 nm: D = S0/4·(λ − λ0⁴/λ³) = 0.023 × (1550 − 1310⁴/1550³) = 0.023 × (1550 − 790.84) = 17.46 ps/(nm·km);
    // Rayleigh 0.9/1.55⁴ = 0.9/5.7720 = 0.15592 dB/km; NA = √(1.4504² − 1.4447²) = 0.12846, V = 2π × 4.1 × 0.12846/1.55 = 2.1350,
    // cut-off = 2π × 4.1e-6 × 0.12846/2.405 = 1376.0 nm; DGD = 0.1 × √100 = 1.0 ps; BER 10⁻⁹ ↔ Q = 5.998 (½·erfc(Q/√2)).
    { name: 'Chromatic dispersion, Rayleigh loss, V number, cut-off, PMD delay and Q for a BER of 10⁻⁹',
      source: 'ITU-T G.652 (dispersion fit); Agrawal, Fiber-Optic Communication Systems §2.2, §2.5, §4.6 (hand calculation)',
      inputs: { L: 100, slack: 0, berReq: 1e-9, dcRatio: 90 },
      expect: [{ kpi: 'D', value: 17.46, tol: 0.001 }, { kpi: 'rayleighLoss', value: 0.15592, tol: 0.001 }, { kpi: 'V', value: 2.1350, tol: 0.001 }, { kpi: 'cutoff', value: 1376.0, tol: 0.001 }, { kpi: 'dgd', value: 1.0, tol: 1e-6 }, { kpi: 'Qreq', value: 5.998, tol: 0.001 }] },
    // Repeatered chain: OSNR = 58 + P − NF − Lspan − 10·log N. The module chooses the smallest N meeting 15 dB for 6000 km at 0.17 dB/km
    // (no slack, no splice or ageing allowance): N = 47 spans → span loss 1020/47 = 21.702 dB,
    // OSNR = 58 + 0 − 4.5 − 21.702 − 16.721 = 15.077 dB (N = 46 gives 22.174 dB and 14.70 dB: not enough). 46 repeaters.
    // Power feed: 46 × 40 V + 1 A × 1 Ω/km × 6000 km + 0.1 V/km × 6000 km = 1840 + 6000 + 600 = 8440 V → 4.22 kV per station.
    { name: 'Amplified chain: OSNR = 58 + P − NF − Lspan − 10·log N, repeater count and power-feed voltage',
      source: 'ITU-T G.977 / G.692 Appendix (OSNR of a chain of identical spans); Ohm\'s law (hand calculation)',
      inputs: { mode: 'repeatered', L: 6000, slack: 0, alpha: 0.17, Ptx: 0, spliceLoss: 0, ageing: 0, dcRatio: 99.5, pmd: 0.05, dts: false, tension: 0 },
      expect: [{ kpi: 'nRepeaters', value: 46, abs: 0 }, { kpi: 'osnr', value: 15.077, tol: 0.001 }, { kpi: 'pfeVoltage', value: 4.22, tol: 0.002 }] },
  ],
};
