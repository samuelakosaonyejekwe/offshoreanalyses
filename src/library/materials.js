// Materials library. Pure data (SI units) so it loads in the page, in workers and under Node.
// Values are nominal specified minimums or handbook-typical figures; `src` names the governing document.
// Project-specific certificates always take precedence — users can add their own records in the app.

const steel = (id, name, fy, fu, src, extra = {}) => ({ id, cat: 'steel', name, E: 207e9, nu: 0.3, rho: 7850, alpha: 1.17e-5, fy: fy * 1e6, fu: fu * 1e6, k: 45, cp: 500, src, ...extra });
const conc = (id, fck, cube, Ecm) => ({ id, cat: 'concrete', name: `Concrete C${fck}/${cube}`, fck: fck * 1e6, fcm: (fck + 8) * 1e6, fctm: (fck <= 50 ? 0.3 * fck ** (2 / 3) : 2.12 * Math.log(1 + (fck + 8) / 10)) * 1e6, E: Ecm * 1e9, nu: 0.2, rho: 2400, rhoReinforced: 2500, alpha: 1e-5, k: 1.6, cp: 900, src: 'EN 1992-1-1 Table 3.1' });

export const MATERIALS = [
  // --- line pipe (API 5L / ISO 3183 PSL2 specified minimum yield / tensile) ---
  steel('api5l-b', 'API 5L Grade B (L245)', 245, 415, 'API 5L / ISO 3183'),
  steel('api5l-x42', 'API 5L X42 (L290)', 290, 415, 'API 5L / ISO 3183'),
  steel('api5l-x52', 'API 5L X52 (L360)', 360, 460, 'API 5L / ISO 3183'),
  steel('api5l-x56', 'API 5L X56 (L390)', 390, 490, 'API 5L / ISO 3183'),
  steel('api5l-x60', 'API 5L X60 (L415)', 415, 520, 'API 5L / ISO 3183'),
  steel('api5l-x65', 'API 5L X65 (L450)', 450, 535, 'API 5L / ISO 3183'),
  steel('api5l-x70', 'API 5L X70 (L485)', 485, 570, 'API 5L / ISO 3183'),
  steel('api5l-x80', 'API 5L X80 (L555)', 555, 625, 'API 5L / ISO 3183'),
  // --- structural steels (EN 10025 / EN 10225 nominal, t ≤ 16 mm) ---
  steel('s235', 'S235 structural steel', 235, 360, 'EN 10025-2', { E: 210e9 }),
  steel('s275', 'S275 structural steel', 275, 430, 'EN 10025-2', { E: 210e9 }),
  steel('s355', 'S355 structural steel', 355, 490, 'EN 10025-2 / EN 10225', { E: 210e9 }),
  steel('s420', 'S420 offshore steel', 420, 500, 'EN 10225', { E: 210e9 }),
  steel('s460', 'S460 offshore steel', 460, 540, 'EN 10225', { E: 210e9 }),
  steel('s690', 'S690 high-strength steel', 690, 770, 'EN 10025-6', { E: 210e9 }),
  steel('a36', 'ASTM A36', 250, 400, 'ASTM A36', { E: 200e9 }),
  steel('a572-50', 'ASTM A572 Grade 50', 345, 450, 'ASTM A572', { E: 200e9 }),
  steel('api2w-50', 'API 2W Grade 50', 345, 448, 'API Spec 2W', { E: 200e9 }),
  steel('api2w-60', 'API 2W Grade 60', 414, 517, 'API Spec 2W', { E: 200e9 }),
  steel('nv-36', 'VL / NV 36 hull steel', 355, 490, 'DNV-OS-B101', { E: 206e9 }),
  steel('ah36', 'AH36 hull steel', 355, 490, 'IACS UR W11', { E: 206e9 }),
  // --- corrosion-resistant alloys ---
  steel('sdss-25cr', 'Super duplex 25Cr (UNS S32750)', 550, 800, 'ASTM A790 / ISO 13628-5', { E: 200e9, rho: 7800, alpha: 1.3e-5, k: 15 }),
  steel('dss-22cr', 'Duplex 22Cr (UNS S31803)', 450, 620, 'ASTM A790', { E: 200e9, rho: 7800, alpha: 1.3e-5, k: 15 }),
  steel('ss316l', 'Stainless 316L', 170, 485, 'ASTM A240', { E: 193e9, rho: 8000, alpha: 1.6e-5, k: 15 }),
  steel('reinf-b500', 'Reinforcing steel B500B', 500, 540, 'EN 10080 / EN 1992-1-1', { E: 200e9 }),
  steel('prestress-1860', 'Prestressing strand Y1860', 1600, 1860, 'EN 10138', { E: 195e9 }),
  steel('armour-wire-hs', 'Flexible-pipe tensile armour wire (high strength, generic)', 1100, 1400, 'generic — replace with manufacturer data', { E: 207e9 }),
  // --- concrete ---
  conc('c25', 25, 30, 31), conc('c30', 30, 37, 33), conc('c35', 35, 45, 34), conc('c40', 40, 50, 35), conc('c45', 45, 55, 36), conc('c50', 50, 60, 37), conc('c60', 60, 75, 39),
  { id: 'cwc-3040', cat: 'concrete', name: 'Concrete weight coating 3040 kg/m³', rho: 3040, fck: 40e6, E: 30e9, nu: 0.2, k: 2.0, cp: 880, waterAbsorption: 0.03, src: 'typical iron-ore aggregate CWC' },
  { id: 'cwc-2400', cat: 'concrete', name: 'Concrete weight coating 2400 kg/m³', rho: 2400, fck: 40e6, E: 28e9, nu: 0.2, k: 1.6, cp: 900, waterAbsorption: 0.04, src: 'typical CWC' },
  { id: 'grout-hs', cat: 'concrete', name: 'High-strength grout (120 MPa class)', rho: 2400, fck: 110e6, E: 50e9, nu: 0.19, src: 'typical offshore grout — replace with supplier data' },
  // --- mooring chain grades: MBL = c · d² (44 − 0.08 d) kN, d in mm (IACS UR W22 / DNV-OS-E302) ---
  ...[['r3', 'R3', 0.0223, 410, 690], ['r3s', 'R3S', 0.0249, 490, 770], ['r4', 'R4', 0.0274, 580, 860], ['r4s', 'R4S', 0.0304, 700, 960], ['r5', 'R5', 0.0320, 760, 1000]]
    .map(([id, g, c, fy, fu]) => ({ id: `chain-${id}`, cat: 'chain', name: `Offshore mooring chain grade ${g}`, grade: g, mblCoeff: c, fy: fy * 1e6, fu: fu * 1e6, E: 207e9, rho: 7850, massCoeffStudless: 19.9e3, massCoeffStud: 21.9e3, eaCoeffStudless: 0.854e11, eaCoeffStud: 1.01e11, src: 'IACS UR W22 / DNV-OS-E302; mass and EA coefficients are indicative (× d² with d in m)' })),
  // --- wire and fibre ropes (indicative; MBL ≈ coeff · d², d in m) ---
  { id: 'wire-spiral', cat: 'wire', name: 'Spiral strand wire rope', mblCoeff: 0.90e9, massCoeff: 4.9e3, eaCoeff: 9.0e10, src: 'indicative — replace with manufacturer catalogue' },
  { id: 'wire-sixstrand', cat: 'wire', name: 'Six-strand wire rope (IWRC)', mblCoeff: 0.63e9, massCoeff: 4.0e3, eaCoeff: 4.5e10, src: 'indicative — replace with manufacturer catalogue' },
  { id: 'rope-polyester', cat: 'rope', name: 'Polyester mooring rope', mblCoeff: 0.17e9, massCoeff: 0.80e3, eaOverMblStatic: 13, eaOverMblDynamic: 27, rhoFibre: 1380, src: 'indicative — replace with manufacturer data' },
  { id: 'rope-hmpe', cat: 'rope', name: 'HMPE mooring rope', mblCoeff: 0.57e9, massCoeff: 0.62e3, eaOverMblStatic: 45, eaOverMblDynamic: 60, rhoFibre: 975, src: 'indicative — replace with manufacturer data' },
  { id: 'rope-nylon', cat: 'rope', name: 'Nylon (polyamide) rope', mblCoeff: 0.14e9, massCoeff: 0.65e3, eaOverMblStatic: 5, eaOverMblDynamic: 9, rhoFibre: 1140, src: 'indicative — replace with manufacturer data' },
  // --- polymers and coatings ---
  { id: 'hdpe', cat: 'polymer', name: 'HDPE (outer sheath)', E: 0.8e9, nu: 0.42, rho: 950, k: 0.41, cp: 1900, fy: 22e6, strainLimit: 0.077, src: 'API 17J typical' },
  { id: 'pa11', cat: 'polymer', name: 'Polyamide PA11 (pressure sheath)', E: 0.35e9, nu: 0.4, rho: 1050, k: 0.33, cp: 2100, fy: 30e6, strainLimit: 0.077, tMax: 90, src: 'API 17J / 17TR2 typical' },
  { id: 'pa12', cat: 'polymer', name: 'Polyamide PA12 (pressure sheath)', E: 0.4e9, nu: 0.4, rho: 1020, k: 0.25, cp: 2100, fy: 35e6, strainLimit: 0.077, tMax: 90, src: 'API 17J typical' },
  { id: 'pvdf', cat: 'polymer', name: 'PVDF (pressure sheath)', E: 0.9e9, nu: 0.38, rho: 1780, k: 0.19, cp: 1200, fy: 40e6, strainLimit: 0.07, tMax: 130, src: 'API 17J typical' },
  { id: 'xlpe', cat: 'insulation', name: 'XLPE insulation', E: 0.2e9, rho: 920, k: 0.286, cp: 2400, epsR: 2.5, tanDelta: 0.001, thetaMax: 90, thetaSc: 250, rhoT: 3.5, src: 'IEC 60287-1-1 / -2-1' },
  { id: 'epr', cat: 'insulation', name: 'EPR insulation', rho: 1300, k: 0.2, epsR: 3.0, tanDelta: 0.005, thetaMax: 90, thetaSc: 250, rhoT: 5.0, src: 'IEC 60287' },
  { id: 'mi-paper', cat: 'insulation', name: 'Mass-impregnated paper (HVDC)', rho: 1100, epsR: 4.0, tanDelta: 0.01, thetaMax: 55, rhoT: 6.0, src: 'IEC 60287' },
  { id: 'pu-stiffener', cat: 'polymer', name: 'Polyurethane (bend stiffener)', E: 0.07e9, nu: 0.45, rho: 1150, strainLimit: 0.1, src: 'generic — replace with supplier data' },
  { id: 'pp-3lpp', cat: 'polymer', name: '3-layer polypropylene coating', rho: 900, k: 0.22, cp: 2000, src: 'typical' },
  { id: 'pp-foam', cat: 'polymer', name: 'PP foam / syntactic insulation', rho: 700, k: 0.17, cp: 2000, src: 'typical wet insulation' },
  { id: 'fbe', cat: 'polymer', name: 'Fusion-bonded epoxy', rho: 1400, k: 0.3, cp: 1500, src: 'typical' },
  { id: 'syntactic-foam', cat: 'polymer', name: 'Syntactic foam buoyancy (2000 m rated)', rho: 480, k: 0.12, src: 'generic — density rises with depth rating' },
  { id: 'gfrp-ud', cat: 'polymer', name: 'Glass/epoxy UD laminate (generic)', E: 40e9, E2: 10e9, G12: 4e9, nu: 0.28, rho: 1900, Xt: 900e6, Xc: 600e6, src: 'generic — replace with tested laminate data' },
  { id: 'cfrp-ud', cat: 'polymer', name: 'Carbon/epoxy UD laminate (generic)', E: 130e9, E2: 9e9, G12: 5e9, nu: 0.3, rho: 1550, Xt: 1800e6, Xc: 1100e6, src: 'generic — replace with tested laminate data' },
  { id: 'rubber-fender', cat: 'polymer', name: 'Fender rubber compound (generic)', E: 4e6, nu: 0.49, rho: 1200, src: 'generic — performance comes from certified fender curves' },
  // --- conductors and metallic cable components ---
  { id: 'cu', cat: 'conductor', name: 'Copper conductor', rho20: 1.7241e-8, alpha20: 3.93e-3, rho: 8890, E: 110e9, cv: 3.45e6, K: 226, beta: 234.5, src: 'IEC 60287-1-1 / IEC 60949' },
  { id: 'al', cat: 'conductor', name: 'Aluminium conductor', rho20: 2.8264e-8, alpha20: 4.03e-3, rho: 2703, E: 69e9, cv: 2.5e6, K: 148, beta: 228, src: 'IEC 60287-1-1 / IEC 60949' },
  { id: 'lead-sheath', cat: 'conductor', name: 'Lead alloy sheath', rho20: 21.4e-8, alpha20: 4.0e-3, rho: 11340, cv: 1.45e6, K: 41, beta: 230, src: 'IEC 60287-1-1 / IEC 60949' },
  { id: 'steel-armour', cat: 'conductor', name: 'Galvanised steel armour wire', rho20: 13.8e-8, alpha20: 4.5e-3, rho: 7850, E: 200e9, cv: 3.8e6, K: 78, beta: 202, fu: 420e6, src: 'IEC 60287-1-1 / IEC 60949' },
  // --- soils (typical ranges; always replace with site investigation) ---
  { id: 'clay-vsoft', cat: 'soil', name: 'Very soft clay', type: 'clay', su: 6e3, suGrad: 1.2e3, gammaSub: 5.0e3, eps50: 0.02, src: 'typical — site data required' },
  { id: 'clay-soft', cat: 'soil', name: 'Soft clay', type: 'clay', su: 20e3, suGrad: 1.5e3, gammaSub: 6.5e3, eps50: 0.02, src: 'typical — site data required' },
  { id: 'clay-firm', cat: 'soil', name: 'Firm clay', type: 'clay', su: 60e3, suGrad: 2.0e3, gammaSub: 8.0e3, eps50: 0.01, src: 'typical — site data required' },
  { id: 'clay-stiff', cat: 'soil', name: 'Stiff clay', type: 'clay', su: 120e3, suGrad: 2.0e3, gammaSub: 9.0e3, eps50: 0.007, src: 'typical — site data required' },
  { id: 'clay-hard', cat: 'soil', name: 'Hard clay', type: 'clay', su: 300e3, suGrad: 0, gammaSub: 10.0e3, eps50: 0.004, src: 'typical — site data required' },
  { id: 'sand-loose', cat: 'soil', name: 'Loose sand', type: 'sand', phiDeg: 30, gammaSub: 8.5e3, Dr: 0.3, d50: 0.2e-3, src: 'typical — site data required' },
  { id: 'sand-medium', cat: 'soil', name: 'Medium-dense sand', type: 'sand', phiDeg: 35, gammaSub: 9.5e3, Dr: 0.55, d50: 0.25e-3, src: 'typical — site data required' },
  { id: 'sand-dense', cat: 'soil', name: 'Dense sand', type: 'sand', phiDeg: 40, gammaSub: 10.0e3, Dr: 0.8, d50: 0.3e-3, src: 'typical — site data required' },
  { id: 'sand-vdense', cat: 'soil', name: 'Very dense sand', type: 'sand', phiDeg: 42, gammaSub: 10.5e3, Dr: 0.95, d50: 0.35e-3, src: 'typical — site data required' },
  { id: 'silt', cat: 'soil', name: 'Silt', type: 'silt', phiDeg: 28, gammaSub: 8.0e3, d50: 0.03e-3, src: 'typical — site data required' },
  // --- rock and fill ---
  { id: 'rock-granite', cat: 'rock', name: 'Armour rock — granite', rho: 2650, src: 'CIRIA Rock Manual typical' },
  { id: 'rock-basalt', cat: 'rock', name: 'Armour rock — basalt', rho: 2900, src: 'CIRIA Rock Manual typical' },
  { id: 'rock-limestone', cat: 'rock', name: 'Armour rock — limestone', rho: 2500, src: 'CIRIA Rock Manual typical' },
  { id: 'concrete-unit', cat: 'rock', name: 'Concrete armour unit', rho: 2350, src: 'typical' },
  // --- fluids ---
  { id: 'seawater', cat: 'fluid', name: 'Seawater (35 PSU, 10 °C)', rho: 1025, mu: 1.4e-3, cp: 3990, k: 0.59, src: 'ITTC 7.5-02-01-03' },
  { id: 'freshwater', cat: 'fluid', name: 'Fresh water (15 °C)', rho: 999, mu: 1.14e-3, cp: 4186, k: 0.59, src: 'handbook' },
  { id: 'air', cat: 'fluid', name: 'Air (15 °C, sea level)', rho: 1.225, mu: 1.81e-5, cp: 1005, k: 0.025, src: 'ISO 2533' },
  { id: 'crude-light', cat: 'fluid', name: 'Light crude oil (35 °API)', rho: 850, mu: 5e-3, cp: 2000, k: 0.13, src: 'typical' },
  { id: 'crude-heavy', cat: 'fluid', name: 'Heavy crude oil (18 °API)', rho: 946, mu: 0.3, cp: 1900, k: 0.13, src: 'typical' },
  { id: 'natgas', cat: 'fluid', name: 'Natural gas (SG 0.65, at 100 bar, 20 °C)', rho: 85, mu: 1.4e-5, cp: 2700, k: 0.045, src: 'typical' },
  { id: 'meg', cat: 'fluid', name: 'Mono-ethylene glycol', rho: 1113, mu: 0.021, cp: 2400, k: 0.25, M: 62.07, src: 'handbook' },
  { id: 'methanol', cat: 'fluid', name: 'Methanol', rho: 792, mu: 5.9e-4, cp: 2500, k: 0.2, M: 32.04, src: 'handbook' },
  { id: 'mgo', cat: 'fluid', name: 'Marine gas oil', rho: 860, lhv: 42.7e6, co2PerKg: 3.206, src: 'IMO MEPC.308(73)' },
];

const byId = new Map(MATERIALS.map((m) => [m.id, m]));
let custom = [];
/** Replace the user-defined records (kept in the project). */
export function setCustomMaterials(list) { custom = Array.isArray(list) ? list : []; }
export const material = (id) => custom.find((m) => m.id === id) || byId.get(id) || null;
export const materialsIn = (cat) => [...MATERIALS, ...custom].filter((m) => !cat || m.cat === cat);
export const CATEGORIES = ['steel', 'concrete', 'chain', 'wire', 'rope', 'polymer', 'insulation', 'conductor', 'soil', 'rock', 'fluid'];
