import { AtomNode, BondEdge, LCAOContribution } from './moleculesData';
import {
  LowdinEHTResult,
  extractLCAOFromLowdin,
  runLowdinGeneralizedEHT,
} from './quantumEngineV3';
import { detectExactPointGroup3D } from './symmetryV3';

export type WalshCoordinateType =
  | 'bend'
  | 'pyramidalize'
  | 'torsion'
  | 'flatten'
  | 'stretch';

export interface WalshCoordinateDescriptor {
  type: WalshCoordinateType;
  label: string;
  unit: string;
  min: number;
  max: number;
  equilibrium: number;
  description: string;
}

export interface WalshFrame {
  value: number;
  atoms: AtomNode[];
  pointGroupId: string;
  schoenflies: string;
  energies: number[]; // All MO energies at this geometry
  homoIndex: number;
  lumoIndex: number;
  totalElectronicEnergy: number;
  dipoleDebye: number;
  frontierLcaoByMoIndex: Record<number, LCAOContribution[]>;
}

export interface WalshScanResult {
  coordinate: WalshCoordinateDescriptor;
  availableCoordinates: WalshCoordinateDescriptor[];
  frames: WalshFrame[];
  moLabels: string[];
  trackedMoIndices: number[]; // Indices of valence/frontier MOs plotted on the correlation diagram
}

function findCentralAtomAndNeighbors(
  atoms: AtomNode[],
  bonds: BondEdge[]
): { centerIdx: number; neighborIndices: number[] } | null {
  if (atoms.length < 3) return null;
  const idToIdx = new Map(atoms.map((a, i) => [a.id, i]));
  const adj: number[][] = atoms.map(() => []);
  for (const b of bonds) {
    const i = idToIdx.get(b.source);
    const j = idToIdx.get(b.target);
    if (i !== undefined && j !== undefined) {
      adj[i].push(j);
      adj[j].push(i);
    }
  }
  let best = -1;
  let bestDeg = 0;
  for (let i = 0; i < atoms.length; i++) {
    if (adj[i].length > bestDeg) {
      bestDeg = adj[i].length;
      best = i;
    }
  }
  if (best < 0 || bestDeg < 2) return null;
  return { centerIdx: best, neighborIndices: adj[best] };
}

/**
 * Finds the best central bond for a torsion scan (e.g., C=C in ethylene, C-C in ethane, O-O in H2O2).
 */
function findTorsionBond(
  atoms: AtomNode[],
  bonds: BondEdge[]
): { aIdx: number; bIdx: number; sideA: number[]; sideB: number[] } | null {
  if (atoms.length < 4) return null;
  const idToIdx = new Map(atoms.map((a, i) => [a.id, i]));
  const adj: number[][] = atoms.map(() => []);
  for (const b of bonds) {
    const i = idToIdx.get(b.source);
    const j = idToIdx.get(b.target);
    if (i !== undefined && j !== undefined) {
      adj[i].push(j);
      adj[j].push(i);
    }
  }
  for (const b of bonds) {
    const i = idToIdx.get(b.source);
    const j = idToIdx.get(b.target);
    if (i === undefined || j === undefined) continue;
    if (adj[i].length >= 2 && adj[j].length >= 2) {
      // BFS on side B without crossing edge (i, j) to check if acyclic torsion
      const visitedB = new Set<number>([j]);
      const queue = [j];
      let cyclic = false;
      while (queue.length > 0) {
        const cur = queue.shift()!;
        for (const nxt of adj[cur]) {
          if (cur === j && nxt === i) continue;
          if (nxt === i) {
            cyclic = true;
            break;
          }
          if (!visitedB.has(nxt)) {
            visitedB.add(nxt);
            queue.push(nxt);
          }
        }
        if (cyclic) break;
      }
      if (!cyclic && visitedB.size >= 2 && visitedB.size < atoms.length - 1) {
        const sideB = Array.from(visitedB);
        const sideA = atoms.map((_, idx) => idx).filter((idx) => !visitedB.has(idx));
        return { aIdx: i, bIdx: j, sideA, sideB };
      }
    }
  }
  return null;
}

/**
 * Auto-detects all chemically valid Walsh perturbation coordinates for the active molecule,
 * ordered with the textbook primary coordinate first.
 */
export function detectAvailableWalshCoordinates(
  atoms: AtomNode[],
  bonds: BondEdge[]
): WalshCoordinateDescriptor[] {
  const available: WalshCoordinateDescriptor[] = [];
  const central = findCentralAtomAndNeighbors(atoms, bonds);

  // 1. Triatomic AB2 bending (H2O, CO2, SO2, NO2, BeH2, XeF2)
  if (atoms.length === 3 && central && central.neighborIndices.length === 2) {
    const c = atoms[central.centerIdx];
    const n1 = atoms[central.neighborIndices[0]];
    const n2 = atoms[central.neighborIndices[1]];
    const v1 = [n1.x3d - c.x3d, n1.y3d - c.y3d, n1.z3d - c.z3d];
    const v2 = [n2.x3d - c.x3d, n2.y3d - c.y3d, n2.z3d - c.z3d];
    const dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    const cosA = Math.max(-1, Math.min(1, dot / (Math.hypot(...v1) * Math.hypot(...v2) || 1)));
    const eqDeg = Math.round((Math.acos(cosA) * 180) / Math.PI);
    available.push({
      type: 'bend',
      label: 'Bond Angle Bend (θ)',
      unit: '°',
      min: 90,
      max: 180,
      equilibrium: Math.max(90, Math.min(180, eqDeg)),
      description: 'Classic AB₂ Walsh correlation from bent C₂v (90°) to linear D∞h (180°).',
    });
  }

  // 2. Tetratomic AX3 umbrella pyramidalization (NH3, BF3, PH3, CH3)
  if (atoms.length === 4 && central && central.neighborIndices.length === 3) {
    const c = atoms[central.centerIdx];
    const ns = central.neighborIndices.map((idx) => atoms[idx]);
    const avgX = (ns[0].x3d + ns[1].x3d + ns[2].x3d) / 3;
    const avgY = (ns[0].y3d + ns[1].y3d + ns[2].y3d) / 3;
    const avgZ = (ns[0].z3d + ns[1].z3d + ns[2].z3d) / 3;
    const hOut = Math.hypot(c.x3d - avgX, c.y3d - avgY, c.z3d - avgZ);
    const rBond = Math.hypot(ns[0].x3d - c.x3d, ns[0].y3d - c.y3d, ns[0].z3d - c.z3d) || 1;
    const eqDeg = Math.round((Math.asin(Math.min(0.8, hOut / rBond)) * 180) / Math.PI);
    available.push({
      type: 'pyramidalize',
      label: 'Umbrella Pyramidalization (φ)',
      unit: '°',
      min: 0,
      max: 34,
      equilibrium: Math.max(0, Math.min(34, eqDeg)),
      description: 'Planar D₃h (0°) to trigonal pyramidal C₃v (34°) umbrella Walsh coordinate.',
    });
  }

  // 3. 4-coordinate ML4 / AX4 square-planar ↔ tetrahedral flattening (CH4, PtCl4, TiCl4, XeF4)
  if (atoms.length === 5 && central && central.neighborIndices.length === 4) {
    // Measure current out-of-plane angle relative to best plane
    const c = atoms[central.centerIdx];
    const ns = central.neighborIndices.map((idx) => atoms[idx]);
    const zSpread = Math.max(...ns.map((a) => Math.abs(a.z3d - c.z3d)));
    const rBond = Math.hypot(ns[0].x3d - c.x3d, ns[0].y3d - c.y3d, ns[0].z3d - c.z3d) || 1;
    const eqDeg = Math.round((Math.asin(Math.min(0.95, zSpread / rBond)) * 180) / Math.PI);
    available.push({
      type: 'flatten',
      label: 'D₄h ↔ Td Distortion (β)',
      unit: '°',
      min: 0,
      max: 35,
      equilibrium: Math.max(0, Math.min(35, eqDeg)),
      description: 'Continuous distortion between square-planar D₄h (0°) and tetrahedral Td (35.3°).',
    });
  }

  // 4. Acyclic torsion twist (C2H4, C2H6, H2O2, etc.)
  if (findTorsionBond(atoms, bonds)) {
    available.push({
      type: 'torsion',
      label: 'Central Bond Torsion (τ)',
      unit: '°',
      min: 0,
      max: 90,
      equilibrium: 0,
      description: 'Internal dihedral twist around the central bond axis (0° equilibrium → 90° twisted).',
    });
  }

  // 5. Universal symmetric bond stretch
  available.push({
    type: 'stretch',
    label: 'Symmetric Stretch (R / R₀)',
    unit: '×',
    min: 0.78,
    max: 1.45,
    equilibrium: 1.0,
    description: 'Uniform symmetric breathing distortion scaling all bond lengths relative to equilibrium.',
  });

  return available;
}

/**
 * Generates a distorted molecular geometry at coordinate `value` for the given Walsh mode.
 */
export function buildPerturbedGeometry(
  atoms: AtomNode[],
  bonds: BondEdge[],
  coordType: WalshCoordinateType,
  value: number
): AtomNode[] {
  const clone = atoms.map((a) => ({ ...a }));

  if (coordType === 'stretch') {
    let cx = 0, cy = 0, cz = 0;
    for (const a of clone) {
      cx += a.x3d;
      cy += a.y3d;
      cz += a.z3d;
    }
    cx /= clone.length || 1;
    cy /= clone.length || 1;
    cz /= clone.length || 1;
    return clone.map((a) => ({
      ...a,
      x3d: cx + (a.x3d - cx) * value,
      y3d: cy + (a.y3d - cy) * value,
      z3d: cz + (a.z3d - cz) * value,
    }));
  }

  const central = findCentralAtomAndNeighbors(clone, bonds);

  if (coordType === 'bend' && clone.length === 3 && central && central.neighborIndices.length === 2) {
    const c = clone[central.centerIdx];
    const n1 = clone[central.neighborIndices[0]];
    const n2 = clone[central.neighborIndices[1]];
    const r1 = Math.hypot(n1.x3d - c.x3d, n1.y3d - c.y3d, n1.z3d - c.z3d) || 1.0;
    const r2 = Math.hypot(n2.x3d - c.x3d, n2.y3d - c.y3d, n2.z3d - c.z3d) || 1.0;
    const halfTheta = ((value * Math.PI) / 180) * 0.5;
    // Place apex at origin pointing +y, ligands symmetric in xy plane
    c.x3d = 0;
    c.y3d = 0.25 * Math.cos(halfTheta);
    c.z3d = 0;
    n1.x3d = -r1 * Math.sin(halfTheta);
    n1.y3d = c.y3d - r1 * Math.cos(halfTheta);
    n1.z3d = 0;
    n2.x3d = r2 * Math.sin(halfTheta);
    n2.y3d = c.y3d - r2 * Math.cos(halfTheta);
    n2.z3d = 0;
    return clone;
  }

  if (coordType === 'pyramidalize' && clone.length === 4 && central && central.neighborIndices.length === 3) {
    const c = clone[central.centerIdx];
    const ns = central.neighborIndices.map((idx) => clone[idx]);
    const rBond =
      ns.reduce(
        (s, a) => s + Math.hypot(a.x3d - c.x3d, a.y3d - c.y3d, a.z3d - c.z3d),
        0
      ) / 3 || 1.05;
    const phi = (value * Math.PI) / 180;
    const rPlane = rBond * Math.cos(phi);
    const yDrop = rBond * Math.sin(phi);
    c.x3d = 0;
    c.y3d = yDrop * 0.5;
    c.z3d = 0;
    for (let k = 0; k < 3; k++) {
      const ang = (k * 2 * Math.PI) / 3;
      ns[k].x3d = rPlane * Math.cos(ang);
      ns[k].y3d = c.y3d - yDrop;
      ns[k].z3d = rPlane * Math.sin(ang);
    }
    return clone;
  }

  if (coordType === 'flatten' && clone.length === 5 && central && central.neighborIndices.length === 4) {
    const c = clone[central.centerIdx];
    const ns = central.neighborIndices.map((idx) => clone[idx]);
    const rBond =
      ns.reduce(
        (s, a) => s + Math.hypot(a.x3d - c.x3d, a.y3d - c.y3d, a.z3d - c.z3d),
        0
      ) / 4 || 1.1;
    const beta = (value * Math.PI) / 180; // 0° = D4h square planar, 35.264° = ideal Td
    const rEq = rBond * Math.cos(beta);
    const zOff = rBond * Math.sin(beta);
    c.x3d = 0;
    c.y3d = 0;
    c.z3d = 0;
    const dirs: [number, number, number][] = [
      [1, 0, 1],
      [-1, 0, 1],
      [0, 1, -1],
      [0, -1, -1],
    ];
    for (let k = 0; k < 4; k++) {
      ns[k].x3d = dirs[k][0] * rEq;
      ns[k].y3d = dirs[k][1] * rEq;
      ns[k].z3d = dirs[k][2] * zOff;
    }
    return clone;
  }

  if (coordType === 'torsion') {
    const tb = findTorsionBond(clone, bonds);
    if (tb) {
      const pA = clone[tb.aIdx];
      const pB = clone[tb.bIdx];
      const ax = pB.x3d - pA.x3d;
      const ay = pB.y3d - pA.y3d;
      const az = pB.z3d - pA.z3d;
      const len = Math.hypot(ax, ay, az) || 1;
      const ux = ax / len,
        uy = ay / len,
        uz = az / len;
      const rad = (value * Math.PI) / 180;
      const cosT = Math.cos(rad);
      const sinT = Math.sin(rad);
      const oneMinusC = 1 - cosT;

      for (const idx of tb.sideB) {
        if (idx === tb.bIdx) continue;
        const a = clone[idx];
        const rx = a.x3d - pB.x3d;
        const ry = a.y3d - pB.y3d;
        const rz = a.z3d - pB.z3d;
        const dot = ux * rx + uy * ry + uz * rz;
        // Rodrigues rotation around bond axis
        a.x3d = pB.x3d + rx * cosT + (uy * rz - uz * ry) * sinT + ux * dot * oneMinusC;
        a.y3d = pB.y3d + ry * cosT + (uz * rx - ux * rz) * sinT + uy * dot * oneMinusC;
        a.z3d = pB.z3d + rz * cosT + (ux * ry - uy * rx) * sinT + uz * dot * oneMinusC;
      }
      return clone;
    }
  }

  return clone;
}

/**
 * Runs a 21-point Löwdin generalized EHT Walsh scan along the selected coordinate.
 */
export function runWalshDiagramScan(
  atoms: AtomNode[],
  bonds: BondEdge[],
  charge = 0,
  preferredType?: WalshCoordinateType,
  numSteps = 19
): WalshScanResult {
  const availableCoordinates = detectAvailableWalshCoordinates(atoms, bonds);
  const coordinate =
    (preferredType && availableCoordinates.find((c) => c.type === preferredType)) ||
    availableCoordinates[0];

  const frames: WalshFrame[] = [];
  let refEht: LowdinEHTResult | null = null;

  for (let s = 0; s < numSteps; s++) {
    const frac = s / Math.max(1, numSteps - 1);
    const val = coordinate.min + frac * (coordinate.max - coordinate.min);
    const distAtoms = buildPerturbedGeometry(atoms, bonds, coordinate.type, val);
    const eht = runLowdinGeneralizedEHT(distAtoms, bonds, charge, false);
    const sym = detectExactPointGroup3D(distAtoms);

    if (!refEht) refEht = eht;

    // Store LCAO for frontier window around HOMO/LUMO for instant live 3D morphing
    const lo = Math.max(0, eht.homoIndex - 5);
    const hi = Math.min(eht.energies.length - 1, Math.max(eht.lumoIndex + 3, eht.homoIndex + 4));
    const frontierLcaoByMoIndex: Record<number, LCAOContribution[]> = {};
    for (let k = lo; k <= hi; k++) {
      frontierLcaoByMoIndex[k] = extractLCAOFromLowdin(eht, k);
    }

    frames.push({
      value: Math.round(val * 1000) / 1000,
      atoms: distAtoms,
      pointGroupId: sym.pointGroupId,
      schoenflies: sym.pointGroup.schoenflies,
      energies: eht.energies,
      homoIndex: eht.homoIndex,
      lumoIndex: eht.lumoIndex,
      totalElectronicEnergy: eht.totalElectronicEnergy,
      dipoleDebye: eht.dipoleMagnitudeDebye,
      frontierLcaoByMoIndex,
    });
  }

  const nMOs = refEht ? refEht.energies.length : 0;
  const homoRef = refEht ? Math.max(0, refEht.homoIndex) : 0;
  // Track up to 10 most chemically relevant valence/frontier MOs
  const startMo = Math.max(0, homoRef - 6);
  const endMo = Math.min(nMOs - 1, homoRef + 3);
  const trackedMoIndices: number[] = [];
  for (let k = startMo; k <= endMo; k++) trackedMoIndices.push(k);

  const moLabels = Array.from({ length: nMOs }, (_, k) => {
    if (k === homoRef) return `MO ${k + 1} (HOMO)`;
    if (k === homoRef + 1) return `MO ${k + 1} (LUMO)`;
    return `MO ${k + 1}`;
  });

  return {
    coordinate,
    availableCoordinates,
    frames,
    moLabels,
    trackedMoIndices,
  };
}
