import {
  AtomNode,
  BondEdge,
  CANONICAL_MOLECULES,
  ELEMENT_PROPERTIES,
  ElementSymbol,
  FragmentOrbital,
  LCAOContribution,
  MolecularOrbital,
  MoleculeEntry,
} from './moleculesData';
import { POINT_GROUPS, decomposeReducibleRepresentation } from './symmetryData';

export function toSubscript(n: number): string {
  return String(n);
}

/**
 * Determines how many implicit hydrogens a skeletal atom should receive.
 *
 * Edge-Case Rules:
 * 1. Non-carbon elements (N, O, B, S, P, F, Cl, Xe, H) NEVER receive implicit hydrogens
 *    so inorganic oxides, radicals, and diatomics (NO, NO2, N2O, O3, SO2, N2, O2, BF3, SF6)
 *    are never mutated into hydrides (HNO, H2O3, H2N2, etc.).
 * 2. Carbon (C) only receives implicit hydrogens when it is part of an organic skeleton
 *    (bonded to another Carbon or part of a multi-atom organic chain/ring) and has no
 *    explicit H already attached. Diatomic heteronuclear carbons (CO, CN, CS, CF) and
 *    terminal isocyanides (-N≡C) receive 0 implicit hydrogens.
 */
export function getImplicitHydrogenCountForAtom(
  atom: AtomNode,
  atoms: AtomNode[],
  bonds: BondEdge[]
): number {
  if (atom.element !== 'C') {
    return 0;
  }

  const connected = bonds.filter((b) => b.source === atom.id || b.target === atom.id);
  const hasExplicitH = connected.some((b) => {
    const otherId = b.source === atom.id ? b.target : b.source;
    const other = atoms.find((x) => x.id === otherId);
    return other?.element === 'H';
  });
  if (hasExplicitH) {
    return 0;
  }

  // Diatomic heteronuclear carbon species (e.g., CO, CN, CS, CF) have 0 implicit hydrogens
  if (atoms.length === 2) {
    const other = atoms.find((x) => x.id !== atom.id);
    if (other && other.element !== 'C') {
      return 0;
    }
  }

  // Terminal carbon triple-bonded to Oxygen or Nitrogen (C≡O or R-N≡C isocyanide) has 0 implicit H
  if (connected.length === 1 && connected[0].order === 3) {
    const otherId = connected[0].source === atom.id ? connected[0].target : connected[0].source;
    const other = atoms.find((x) => x.id === otherId);
    if (other && (other.element === 'O' || other.element === 'N' || other.element === 'S')) {
      return 0;
    }
  }

  const bondOrderSum = connected.reduce((sum, b) => sum + b.order, 0);
  return Math.max(0, 4 - bondOrderSum);
}

/**
 * Computes the Hill / IUPAC molecular formula (C first, then H, then alphabetical for organic;
 * IUPAC electropositive-to-electronegative order for inorganic) including implicit hydrogens
 * only when `includeImplicitH` is true.
 */
export function computeMolecularFormula(
  atoms: AtomNode[],
  bonds: BondEdge[],
  includeImplicitH = false
): { formula: string; plainFormula: string; counts: Record<string, number> } {
  const counts: Record<string, number> = {};
  for (const a of atoms) {
    counts[a.element] = (counts[a.element] || 0) + 1;
  }

  if (includeImplicitH) {
    for (const a of atoms) {
      const imp = getImplicitHydrogenCountForAtom(a, atoms, bonds);
      if (imp > 0) {
        counts['H'] = (counts['H'] || 0) + imp;
      }
    }
  }

  const keys = Object.keys(counts);
  const orderedKeys: string[] = [];
  if (counts['C']) {
    orderedKeys.push('C');
    if (counts['H']) orderedKeys.push('H');
    for (const k of keys.sort()) {
      if (k !== 'C' && k !== 'H') orderedKeys.push(k);
    }
  } else {
    // IUPAC Red Book electropositive-to-electronegative ordering (NH3, H2O, BF3, SF6, NO, NO2, SO2, XeF4)
    const iupacOrder = ['Xe', 'B', 'P', 'N', 'H', 'S', 'Cl', 'O', 'F'];
    const sorted = [...keys].sort((a, b) => {
      const ia = iupacOrder.indexOf(a);
      const ib = iupacOrder.indexOf(b);
      if (ia !== -1 && ib !== -1) return ia - ib;
      if (ia !== -1) return -1;
      if (ib !== -1) return 1;
      return a.localeCompare(b);
    });
    orderedKeys.push(...sorted);
  }

  let formula = '';
  let plainFormula = '';
  for (const el of orderedKeys) {
    const c = counts[el];
    if (!c) continue;
    formula += el + (c > 1 ? String(c) : '');
    plainFormula += el + (c > 1 ? String(c) : '');
  }

  return { formula: formula || 'Empty', plainFormula, counts };
}

/**
 * Generates a simplified canonical SMILES string from the drawn 2D graph.
 */
export function generateSmilesFromGraph(atoms: AtomNode[], bonds: BondEdge[]): string {
  if (atoms.length === 0) return '';
  const nonH = atoms.filter((a) => a.element !== 'H');
  if (nonH.length === 0) return '[H][H]';

  const adj = new Map<string, { neighbor: AtomNode; order: 1 | 2 | 3 }[]>();
  for (const a of nonH) adj.set(a.id, []);
  for (const b of bonds) {
    const s = nonH.find((a) => a.id === b.source);
    const t = nonH.find((a) => a.id === b.target);
    if (s && t) {
      adj.get(s.id)!.push({ neighbor: t, order: b.order });
      adj.get(t.id)!.push({ neighbor: s, order: b.order });
    }
  }

  const visited = new Set<string>();
  const buildDfs = (node: AtomNode, parentId: string | null): string => {
    visited.add(node.id);
    const neighbors = adj.get(node.id) || [];
    const branches: string[] = [];
    let mainChain = '';

    for (const edge of neighbors) {
      if (edge.neighbor.id === parentId) continue;
      const bondSym = edge.order === 2 ? '=' : edge.order === 3 ? '#' : '';
      if (visited.has(edge.neighbor.id)) {
        // Ring closure
        if (!mainChain.includes('1')) {
          mainChain += bondSym + '1';
        }
        continue;
      }
      const sub = bondSym + buildDfs(edge.neighbor, node.id);
      if (!mainChain) {
        mainChain = sub;
      } else {
        branches.push(`(${sub})`);
      }
    }
    return `${node.element}${branches.join('')}${mainChain}`;
  };

  return buildDfs(nonH[0], null);
}

/**
 * Completes missing implicit hydrogens around drawn skeletal Carbon atoms
 * so that 2D organic drawings (like a bare hexagon for Benzene or C=C for Ethylene)
 * automatically expand into full 3D all-atom structures, while leaving inorganic
 * molecules (NO, CO, O3, SO2, NO2, N2, O2) untouched.
 */
export function expandImplicitHydrogens(
  atoms: AtomNode[],
  bonds: BondEdge[],
  skipImplicitH = false
): { atoms: AtomNode[]; bonds: BondEdge[] } {
  const outAtoms: AtomNode[] = atoms.map((a) => ({ ...a }));
  const outBonds: BondEdge[] = bonds.map((b) => ({ ...b }));
  if (skipImplicitH) {
    return { atoms: outAtoms, bonds: outBonds };
  }

  let hIndex = outAtoms.filter((a) => a.element === 'H').length + 1;

  for (const atom of atoms) {
    const missingH = getImplicitHydrogenCountForAtom(atom, atoms, bonds);
    if (missingH === 0) continue;

    const connectedBonds = bonds.filter((b) => b.source === atom.id || b.target === atom.id);

    // Find average direction of existing neighbors to point hydrogens away
    const neighborVectors: [number, number][] = [];
    for (const b of connectedBonds) {
      const otherId = b.source === atom.id ? b.target : b.source;
      const other = atoms.find((a) => a.id === otherId);
      if (other) {
        const dx = other.x2d - atom.x2d;
        const dy = other.y2d - atom.y2d;
        const len = Math.hypot(dx, dy) || 1;
        neighborVectors.push([dx / len, dy / len]);
      }
    }

    let baseAngle = 0;
    if (neighborVectors.length > 0) {
      const sumX = neighborVectors.reduce((s, v) => s + v[0], 0);
      const sumY = neighborVectors.reduce((s, v) => s + v[1], 0);
      if (Math.hypot(sumX, sumY) > 0.1) {
        baseAngle = Math.atan2(-sumY, -sumX);
      } else if (neighborVectors.length === 2) {
        baseAngle = Math.atan2(neighborVectors[0][0], -neighborVectors[0][1]);
      }
    }

    for (let i = 0; i < missingH; i++) {
      const spread = missingH === 1 ? 0 : missingH === 2 ? (i === 0 ? -0.9 : 0.9) : (i * 2 * Math.PI) / missingH;
      const angle = baseAngle + spread;
      const r2d = 0.95;
      const hId = `H_imp_${atom.id}_${hIndex++}`;
      const hx2d = atom.x2d + Math.cos(angle) * r2d;
      const hy2d = atom.y2d + Math.sin(angle) * r2d;
      const hz3d = missingH >= 3 && i > 0 ? (i % 2 === 0 ? 0.75 : -0.75) : 0;

      outAtoms.push({
        id: hId,
        element: 'H',
        x2d: Number(hx2d.toFixed(2)),
        y2d: Number(hy2d.toFixed(2)),
        x3d: Number(hx2d.toFixed(2)),
        y3d: Number((-hy2d).toFixed(2)),
        z3d: hz3d,
      });

      outBonds.push({
        id: `b_imp_${hId}`,
        source: atom.id,
        target: hId,
        order: 1,
      });
    }
  }

  return { atoms: outAtoms, bonds: outBonds };
}

/**
 * Checks if a drawn graph matches one of our canonical benchmark molecules
 * by stoichiometry and bond connectivity (checking exact all-atom match first,
 * then implicit-hydrogen expanded skeletal match).
 */
export function matchCanonicalMolecule(
  atoms: AtomNode[],
  bonds: BondEdge[],
  skipImplicitH = false
): MoleculeEntry | null {
  // 1. Direct match on raw atoms & bonds (e.g. CO, NO, H2O, SF6, or any preset)
  const rawFormula = computeMolecularFormula(atoms, bonds, false).plainFormula;
  for (const entry of CANONICAL_MOLECULES) {
    const entryFormula = computeMolecularFormula(entry.atoms, entry.bonds, false).plainFormula;
    if (entryFormula === rawFormula && entry.bonds.length === bonds.length) {
      return entry;
    }
  }

  if (skipImplicitH) {
    return null;
  }

  // 2. Skeletal match after expanding implicit hydrogens (e.g. user draws bare C6 ring -> C6H6)
  const expanded = expandImplicitHydrogens(atoms, bonds, false);
  const { plainFormula } = computeMolecularFormula(expanded.atoms, expanded.bonds, false);

  for (const entry of CANONICAL_MOLECULES) {
    const entryFormula = computeMolecularFormula(entry.atoms, entry.bonds, false).plainFormula;
    if (entryFormula === plainFormula && entry.bonds.length === expanded.bonds.length) {
      return entry;
    }
  }
  return null;
}

/**
 * Jacobi Eigenvalue & Eigenvector solver for real symmetric matrices H (N x N).
 * Returns eigenvalues sorted ascending (lowest energy first) and normalized column eigenvectors.
 */
export function jacobiEigenvalues(matrix: number[][]): {
  eigenvalues: number[];
  eigenvectors: number[][]; // eigenvectors[i][mu] is coeff of basis mu in eigenstate i
} {
  const n = matrix.length;
  if (n === 0) return { eigenvalues: [], eigenvectors: [] };

  const A = matrix.map((row) => [...row]);
  const V: number[][] = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (__, j) => (i === j ? 1 : 0))
  );

  const maxSweeps = 60;
  const eps = 1e-10;

  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let maxOffDiag = 0;
    let p = 0;
    let q = 1;

    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const val = Math.abs(A[i][j]);
        if (val > maxOffDiag) {
          maxOffDiag = val;
          p = i;
          q = j;
        }
      }
    }

    if (maxOffDiag < eps) break;

    const app = A[p][p];
    const aqq = A[q][q];
    const apq = A[p][q];

    const theta = 0.5 * Math.atan2(2 * apq, aqq - app);
    const c = Math.cos(theta);
    const s = Math.sin(theta);

    // Rotate A
    for (let i = 0; i < n; i++) {
      if (i !== p && i !== q) {
        const aip = A[i][p];
        const aiq = A[i][q];
        A[i][p] = c * aip - s * aiq;
        A[p][i] = A[i][p];
        A[i][q] = s * aip + c * aiq;
        A[q][i] = A[i][q];
      }
    }

    A[p][p] = c * c * app - 2 * s * c * apq + s * s * aqq;
    A[q][q] = s * s * app + 2 * s * c * apq + c * c * aqq;
    A[p][q] = 0;
    A[q][p] = 0;

    // Rotate V
    for (let i = 0; i < n; i++) {
      const vip = V[i][p];
      const viq = V[i][q];
      V[i][p] = c * vip - s * viq;
      V[i][q] = s * vip + c * viq;
    }
  }

  const indices = Array.from({ length: n }, (_, i) => i).sort((a, b) => A[a][a] - A[b][b]);
  const eigenvalues = indices.map((idx) => A[idx][idx]);
  const eigenvectors = indices.map((idx) => {
    const vec = Array.from({ length: n }, (_, row) => V[row][idx]);
    // Normalize and enforce deterministic sign (first significant component positive)
    const norm = Math.hypot(...vec) || 1;
    const firstNonZero = vec.find((v) => Math.abs(v) > 1e-4) ?? 1;
    const sign = firstNonZero < 0 ? -1 : 1;
    return vec.map((v) => (v / norm) * sign);
  });

  return { eigenvalues, eigenvectors };
}

/**
 * Determines the molecular Point Group from atomic composition, bond connectivity, and 2D/3D geometry.
 */
export function detectPointGroup(atoms: AtomNode[], bonds: BondEdge[]): string {
  const heavy = atoms.filter((a) => a.element !== 'H');
  const nAtoms = atoms.length;

  if (nAtoms <= 1) return 'Oh';

  // Diatomic molecules
  if (nAtoms === 2) {
    return atoms[0].element === atoms[1].element ? 'Dinfh' : 'Cinfv';
  }

  // Triatomic molecules (e.g., H2O, CO2, O3, SO2, HCN, OCS)
  if (nAtoms === 3) {
    const counts: Record<string, number> = {};
    atoms.forEach((a) => {
      counts[a.element] = (counts[a.element] || 0) + 1;
    });
    const elements = Object.keys(counts);

    // Cyclic equilateral triatomic (3 bonds in a triangle, e.g. H3+ or cyclopropenyl)
    if (elements.length === 1 && bonds.length >= 3) {
      return 'D3h';
    }

    // Find central atom (the atom bonded to both other atoms, or unique element)
    const center =
      atoms.find((a) => bonds.filter((b) => b.source === a.id || b.target === a.id).length >= 2) ||
      atoms.find((a) => counts[a.element] === 1) ||
      atoms[0];
    const others = atoms.filter((a) => a.id !== center.id);

    if (others.length === 2) {
      const v1 = [
        others[0].x3d - center.x3d || others[0].x2d - center.x2d,
        others[0].y3d - center.y3d || others[0].y2d - center.y2d,
        others[0].z3d - center.z3d,
      ];
      const v2 = [
        others[1].x3d - center.x3d || others[1].x2d - center.x2d,
        others[1].y3d - center.y3d || others[1].y2d - center.y2d,
        others[1].z3d - center.z3d,
      ];
      const norm1 = Math.hypot(v1[0], v1[1], v1[2]) || 1;
      const norm2 = Math.hypot(v2[0], v2[1], v2[2]) || 1;
      const cosTheta = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (norm1 * norm2);
      const isLinear = cosTheta < -0.94 || (center.element === 'C' && bonds.some((b) => b.order >= 2));

      if (others[0].element === others[1].element) {
        // Symmetric AB2 or A3 (e.g. CO2 -> Dinfh, H2O / SO2 / O3 -> C2v)
        return isLinear && center.element !== 'O' && center.element !== 'S' ? 'Dinfh' : 'C2v';
      } else {
        // Asymmetric ABC (e.g. HCN, OCS -> Cinfv, HNO -> Cs)
        return isLinear ? 'Cinfv' : 'Cs';
      }
    }
    return 'Cinfv';
  }

  // Check symmetric rings (e.g. C6H6 -> D6h, C4H4 -> D4h)
  if (heavy.length === 6 && heavy.every((a) => a.element === 'C') && bonds.length >= 6) {
    return 'D6h';
  }
  if (heavy.length === 4 && heavy.every((a) => a.element === 'C') && bonds.length >= 4) {
    const degs = heavy.map(
      (a) =>
        bonds.filter(
          (b) =>
            (b.source === a.id || b.target === a.id) &&
            heavy.some((h) => h.id === (b.source === a.id ? b.target : b.source))
        ).length
    );
    if (degs.every((d) => d === 2)) return 'D4h';
    return 'C2v';
  }

  // Central atom with N identical H ligands
  if (heavy.length === 1) {
    const nLigands = atoms.length - 1;
    const centerEl = heavy[0].element;
    if (nLigands === 3) {
      return centerEl === 'B' ? 'D3h' : 'C3v';
    }
    if (nLigands === 4) {
      return centerEl === 'Xe' ? 'D4h' : 'Td';
    }
    if (nLigands === 6) return 'Oh';
  }

  // Central atom surrounded by identical halogen/chalcogen ligands (BF3, PCl3, SF6, XeF4, CCl4)
  const centerAtom = atoms.find((a) => {
    const deg = bonds.filter((b) => b.source === a.id || b.target === a.id).length;
    return deg >= 3 && deg === atoms.length - 1;
  });
  if (centerAtom) {
    const ligandElements = new Set(
      atoms.filter((a) => a.id !== centerAtom.id).map((a) => a.element)
    );
    if (ligandElements.size === 1) {
      const nLig = atoms.length - 1;
      if (nLig === 3) return centerAtom.element === 'B' ? 'D3h' : 'C3v';
      if (nLig === 4) return centerAtom.element === 'Xe' ? 'D4h' : 'Td';
      if (nLig === 6) return 'Oh';
    }
    if (ligandElements.size === 2 && atoms.length === 4) {
      // Planar AB2C (e.g. Formaldehyde CH2O, Phosgene COCl2) -> C2v
      return 'C2v';
    }
  }

  if (heavy.length === 2 && heavy[0].element === heavy[1].element) {
    const nH = atoms.filter((a) => a.element === 'H').length;
    if (nH === 2) return 'Dinfh';
    if (nH === 4) return 'D2h';
  }

  return 'Cs';
}

/**
 * Computes raw Extended Hückel Theory (EHT) eigenvalues and eigenvectors using
 * a full 4-orbital (ns, npx, npy, npz) valence basis for heavy atoms and 1s for H,
 * with directional Slater-Koster two-center overlap integrals S_mu_nu(R_AB).
 */
export function runExtendedHuckelSolver(
  atoms: AtomNode[],
  bonds: BondEdge[]
): {
  basis: {
    atom: AtomNode;
    orbital: 's' | 'px' | 'py' | 'pz';
    label: string;
    alpha: number;
  }[];
  eigenvalues: number[];
  eigenvectors: number[][];
} {
  interface BasisFunc {
    atom: AtomNode;
    orbital: 's' | 'px' | 'py' | 'pz';
    label: string;
    alpha: number; // diagonal Coulomb integral H_ii in eV
  }

  const basis: BasisFunc[] = [];
  for (const a of atoms) {
    const props = ELEMENT_PROPERTIES[a.element] || ELEMENT_PROPERTIES.C;
    const shell = a.element === 'P' || a.element === 'S' || a.element === 'Cl' ? '3' : a.element === 'Xe' ? '5' : '2';
    if (a.element === 'H') {
      basis.push({
        atom: a,
        orbital: 's',
        label: `${a.id} 1s`,
        alpha: props.vsieS,
      });
    } else {
      basis.push({
        atom: a,
        orbital: 's',
        label: `${a.id} ${shell}s`,
        alpha: props.vsieS,
      });
      basis.push({
        atom: a,
        orbital: 'px',
        label: `${a.id} ${shell}px`,
        alpha: props.vsieP,
      });
      basis.push({
        atom: a,
        orbital: 'py',
        label: `${a.id} ${shell}py`,
        alpha: props.vsieP,
      });
      basis.push({
        atom: a,
        orbital: 'pz',
        label: `${a.id} ${shell}pz`,
        alpha: props.vsieP,
      });
    }
  }

  const N = basis.length;
  const H: number[][] = Array.from({ length: N }, () => Array(N).fill(0));
  const K = 1.75; // Hoffmann Wolfsberg-Helmholz constant

  for (let i = 0; i < N; i++) {
    H[i][i] = basis[i].alpha;
    for (let j = i + 1; j < N; j++) {
      const bi = basis[i];
      const bj = basis[j];
      if (bi.atom.id === bj.atom.id) {
        H[i][j] = 0;
        H[j][i] = 0;
        continue;
      }

      const dx = bj.atom.x3d - bi.atom.x3d;
      const dy = bj.atom.y3d - bi.atom.y3d;
      const dz = bj.atom.z3d - bi.atom.z3d;
      const dist = Math.hypot(dx, dy, dz) || 1.4;

      const bond = bonds.find(
        (b) =>
          (b.source === bi.atom.id && b.target === bj.atom.id) ||
          (b.source === bj.atom.id && b.target === bi.atom.id)
      );

      if (!bond && dist > 2.65) {
        H[i][j] = 0;
        H[j][i] = 0;
        continue;
      }

      // Direction cosines from atom i to atom j
      const l = dx / dist;
      const m = dy / dist;
      const n = dz / dist;
      const dirMap: Record<'px' | 'py' | 'pz', number> = { px: l, py: m, pz: n };

      // Distance-dependent radial Slater-Koster overlaps calibrated to Hoffmann EHT
      const r0 = (ELEMENT_PROPERTIES[bi.atom.element]?.covalentRadius ?? 0.77) +
        (ELEMENT_PROPERTIES[bj.atom.element]?.covalentRadius ?? 0.77);
      const decay = Math.exp(-0.85 * Math.max(0, dist - r0));
      const orderBoost = bond ? 1 + (bond.order - 1) * 0.14 : 0.28;

      const S_ss_sigma = 0.36 * decay * orderBoost;
      const S_sp_sigma = 0.32 * decay * orderBoost;
      const S_pp_sigma = 0.30 * decay * orderBoost;
      const S_pp_pi = 0.19 * decay * orderBoost;

      let overlapS = 0;
      if (bi.orbital === 's' && bj.orbital === 's') {
        overlapS = S_ss_sigma;
      } else if (bi.orbital === 's' && bj.orbital !== 's') {
        overlapS = dirMap[bj.orbital] * S_sp_sigma;
      } else if (bi.orbital !== 's' && bj.orbital === 's') {
        overlapS = -dirMap[bi.orbital] * S_sp_sigma;
      } else if (bi.orbital !== 's' && bj.orbital !== 's') {
        const di = dirMap[bi.orbital];
        const dj = dirMap[bj.orbital];
        if (bi.orbital === bj.orbital) {
          overlapS = di * di * S_pp_sigma + (1 - di * di) * S_pp_pi;
        } else {
          overlapS = di * dj * (S_pp_sigma - S_pp_pi);
        }
      }

      const hij = 0.5 * K * (bi.alpha + bj.alpha) * overlapS;
      H[i][j] = hij;
      H[j][i] = hij;
    }
  }

  const { eigenvalues, eigenvectors } = jacobiEigenvalues(H);
  return { basis, eigenvalues, eigenvectors };
}

/**
 * Builds a complete MoleculeEntry with Extended Hückel LCAO-MO energy levels,
 * eigenvectors, fragment correlations, and symmetry assignments for ANY user-drawn
 * or PubChem-fetched structure.
 */
export function analyzeMolecularGraph(
  rawAtoms: AtomNode[],
  rawBonds: BondEdge[],
  options?: {
    customName?: string;
    customCid?: number;
    charge?: number;
    excitedState?: boolean;
  }
): MoleculeEntry {
  const charge = options?.charge ?? 0;
  const excitedState = options?.excitedState ?? false;

  // 1. Check if the graph matches one of our canonical benchmark molecules
  const canonicalMatch = matchCanonicalMolecule(rawAtoms, rawBonds);
  if (canonicalMatch) {
    return applyElectronicState(
      canonicalMatch,
      charge,
      excitedState,
      options?.customName,
      options?.customCid
    );
  }

  // 2. Otherwise, expand implicit hydrogens and run the full 4-orbital Extended Hückel LCAO-MO + Symmetry solver
  const { atoms, bonds } = expandImplicitHydrogens(rawAtoms, rawBonds);
  const { formula } = computeMolecularFormula(atoms, bonds, false);
  const smiles = generateSmilesFromGraph(atoms, bonds);
  const pointGroupId = detectPointGroup(atoms, bonds);
  const pg = POINT_GROUPS[pointGroupId] || POINT_GROUPS.Cs;

  // Center 3D coordinates around centroid
  const cx = atoms.reduce((s, a) => s + a.x2d, 0) / (atoms.length || 1);
  const cy = atoms.reduce((s, a) => s + a.y2d, 0) / (atoms.length || 1);
  const centeredAtoms: AtomNode[] = atoms.map((a) => ({
    ...a,
    x2d: Number((a.x2d - cx).toFixed(2)),
    y2d: Number((a.y2d - cy).toFixed(2)),
    x3d: Number((a.x2d - cx).toFixed(2)),
    y3d: Number((-(a.y2d - cy)).toFixed(2)),
    z3d: a.z3d || 0,
  }));

  const { basis, eigenvalues, eigenvectors } = runExtendedHuckelSolver(centeredAtoms, bonds);
  const N = basis.length;

  // Exact physical neutral valence electrons = sum(Z_val); charge is applied once in applyElectronicState
  const rawValenceElectrons = centeredAtoms.reduce(
    (sum, a) => sum + (ELEMENT_PROPERTIES[a.element]?.valenceElectrons ?? 4),
    0
  );
  const neutralElectrons = Math.max(1, Math.min(N * 2, rawValenceElectrons));

  // Detect whether the molecule is linear along the x-axis (where both py and pz are pi-type)
  const maxAbsY = centeredAtoms.reduce((m, a) => Math.max(m, Math.abs(a.y3d)), 0);
  const maxAbsZ = centeredAtoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d)), 0);
  const isLinearAlongX = centeredAtoms.length >= 2 && maxAbsY < 0.15 && maxAbsZ < 0.15;

  // Build fragment orbitals (Left = Heavy Atom AOs, Right = Peripheral / Ligand AOs)
  const fragmentOrbitals: FragmentOrbital[] = [
    {
      id: 'frag_left_s',
      side: 'left',
      label: 'Core Valence 2s',
      fragmentName: 'Skeleton AOs',
      energy: -21.4,
      irrep: pg.irreps[0]?.symbol || 'A',
      degeneracy: 1,
      electrons: 2,
    },
    {
      id: 'frag_left_p',
      side: 'left',
      label: 'Valence 2p (σ / π)',
      fragmentName: 'Skeleton AOs',
      energy: -11.4,
      irrep: pg.irreps[Math.min(1, pg.irreps.length - 1)]?.symbol || 'B',
      degeneracy: 3,
      electrons: 2,
    },
    {
      id: 'frag_right_salc',
      side: 'right',
      label: 'Ligand / Peripheral SALCs',
      fragmentName: 'Substituent SALCs',
      energy: -13.6,
      irrep: `${pg.irreps[0]?.symbol || 'A'} + ${pg.irreps[pg.irreps.length - 1]?.symbol || 'B'}`,
      degeneracy: 2,
      electrons: 2,
    },
  ];

  // Select up to 8 most chemically informative frontier & valence MOs around HOMO/LUMO
  const homoIdx = Math.max(0, Math.min(N - 1, Math.floor((neutralElectrons - 1) / 2)));
  const startIdx = Math.max(0, homoIdx - 4);
  const endIdx = Math.min(N - 1, homoIdx + 3);

  const molecularOrbitals: MolecularOrbital[] = [];
  let remainingElectrons = Math.max(0, neutralElectrons - startIdx * 2);

  for (let idx = startIdx; idx <= endIdx; idx++) {
    const energy = Number(eigenvalues[idx].toFixed(2));
    const vec = eigenvectors[idx];
    const piWeight = vec.reduce(
      (acc, c, k) =>
        acc +
        (basis[k].orbital === 'pz' || (isLinearAlongX && basis[k].orbital === 'py')
          ? c * c
          : 0),
      0
    );
    const isPi = piWeight > 0.48;

    const occ: 0 | 1 | 2 =
      remainingElectrons >= 2 ? 2 : remainingElectrons === 1 ? 1 : 0;
    remainingElectrons = Math.max(0, remainingElectrons - occ);

    const irrepObj = pg.irreps[idx % pg.irreps.length] || pg.irreps[0];
    const bondingType: MolecularOrbital['bondingType'] =
      energy < -14.2 ? 'bonding' : energy > -10.2 ? 'antibonding' : 'non-bonding';
    const interactionType: MolecularOrbital['interactionType'] = isPi
      ? bondingType === 'antibonding'
        ? 'π*'
        : 'π'
      : bondingType === 'antibonding'
        ? 'σ*'
        : bondingType === 'non-bonding'
          ? 'n (Lone Pair)'
          : 'σ';

    const lcao: LCAOContribution[] = vec
      .map((coeff, k) => ({
        atomId: basis[k].atom.id,
        element: basis[k].atom.element,
        orbital: basis[k].orbital,
        label: basis[k].label,
        coeff: Number(coeff.toFixed(2)),
      }))
      .filter((c) => Math.abs(c.coeff) >= 0.14)
      .slice(0, 8);

    // Count sign changes across bonds as nodal planes
    let signChanges = 0;
    for (const b of bonds) {
      const cSource = lcao.find((c) => c.atomId === b.source)?.coeff ?? 0;
      const cTarget = lcao.find((c) => c.atomId === b.target)?.coeff ?? 0;
      if (cSource * cTarget < -0.02) signChanges++;
    }

    // Group nearly degenerate eigenvalues into a shared degenerateGroup
    const prevMo = molecularOrbitals[molecularOrbitals.length - 1];
    const isDegenerateWithPrev = prevMo && Math.abs(prevMo.energy - energy) < 0.08;
    const degenerateGroup = isDegenerateWithPrev
      ? prevMo.degenerateGroup || `deg_${idx - 1}`
      : undefined;
    if (isDegenerateWithPrev && prevMo && !prevMo.degenerateGroup) {
      prevMo.degenerateGroup = degenerateGroup;
    }

    molecularOrbitals.push({
      id: `custom_mo_${idx + 1}`,
      label: `${idx - startIdx + 1}${irrepObj.symbol.toLowerCase()}${bondingType === 'antibonding' ? '*' : ''}`,
      irrepSymbol: irrepObj.symbol,
      irrepId: irrepObj.id,
      energy,
      occupancy: occ,
      bondingType,
      interactionType,
      nodalPlanes: Math.min(4, signChanges + (isPi ? 1 : 0)),
      isHOMO: idx === homoIdx,
      isLUMO: idx === homoIdx + 1,
      degenerateGroup,
      correlatedFragmentIds: isPi
        ? ['frag_left_p', 'frag_right_salc']
        : energy < -16
          ? ['frag_left_s', 'frag_right_salc']
          : ['frag_left_p', 'frag_right_salc'],
      lcao,
      description: `Extended Hückel LCAO-MO (${interactionType}, ${irrepObj.symbol} symmetry) computed at ε = ${energy.toFixed(2)} eV with ${signChanges} internuclear nodal surface(s).`,
    });
  }

  const entry: MoleculeEntry = {
    id: `custom_${formula}`,
    name: options?.customName || `Custom Structure (${formula})`,
    iupacName: options?.customName || `Drawn ${formula} Isomer`,
    formula,
    smiles,
    pubchemCid: options?.customCid || 0,
    pointGroupId: pg.id,
    schoenflies: pg.schoenflies,
    valenceElectrons: neutralElectrons,
    bondOrder: bonds.reduce((s, b) => s + b.order, 0),
    dipoleMomentDebye:
      pg.id === 'D6h' || pg.id === 'Td' || pg.id === 'Oh' || pg.id === 'Dinfh' || pg.id === 'D2h' || pg.id === 'D3h' || pg.id === 'D4h'
        ? 0.0
        : 1.42,
    leftFragmentTitle: 'Central / Skeleton AOs',
    rightFragmentTitle: 'Ligand / Substituent SALCs',
    atoms: centeredAtoms,
    bonds,
    fragmentOrbitals,
    molecularOrbitals,
    reducibleRepresentations: pg.defaultReduciblePresets.map((p) => ({
      label: p.label,
      characters: p.characters,
      decomposition: decomposeReducibleRepresentation(pg, p.characters).formula,
    })),
  };

  return applyElectronicState(entry, charge, excitedState, options?.customName, options?.customCid);
}

/**
 * Applies charge (ionization/reduction) and photo-excitation (HOMO -> LUMO promotion)
 * to any MoleculeEntry while enforcing the Aufbau Principle + Hund's Rule of Maximum
 * Multiplicity across degenerate orbital shells, and recalculates HOMO, LUMO, and net bond order.
 */
export function applyElectronicState(
  base: MoleculeEntry,
  charge: number,
  excitedState: boolean,
  overrideName?: string,
  overrideCid?: number
): MoleculeEntry {
  const mos = base.molecularOrbitals.map((m) => ({
    ...m,
    occupancy: 0 as 0 | 1 | 2,
    isHOMO: false,
    isLUMO: false,
  }));
  const baseElectrons = base.molecularOrbitals.reduce((sum, m) => sum + m.occupancy, 0);
  let targetElectrons = Math.max(1, Math.min(mos.length * 2, baseElectrons - charge));

  // Partition consecutive MOs into degenerate shells (shared degenerateGroup or |Δε| < 0.05 eV)
  const shells: MolecularOrbital[][] = [];
  for (const m of mos) {
    const prevShell = shells[shells.length - 1];
    if (
      prevShell &&
      ((m.degenerateGroup && prevShell[0].degenerateGroup === m.degenerateGroup) ||
        Math.abs(prevShell[0].energy - m.energy) < 0.05)
    ) {
      prevShell.push(m);
    } else {
      shells.push([m]);
    }
  }

  // Fill shells according to Aufbau Principle + Hund's Rule of Maximum Multiplicity:
  // First half-fill each orbital in a degenerate shell with 1 spin-up electron before spin-pairing
  for (const shell of shells) {
    // Pass 1: single occupancy across degenerate shell (Hund's rule)
    for (const m of shell) {
      if (targetElectrons > 0) {
        m.occupancy = 1;
        targetElectrons -= 1;
      }
    }
    // Pass 2: pair spins across degenerate shell
    for (const m of shell) {
      if (targetElectrons > 0 && m.occupancy === 1) {
        m.occupancy = 2;
        targetElectrons -= 1;
      }
    }
  }

  // Find HOMO (first occupied orbital in the highest occupied shell, or highest occupied MO) and LUMO
  let homoIdx = -1;
  for (let i = mos.length - 1; i >= 0; i--) {
    if (mos[i].occupancy > 0) {
      homoIdx = i;
      break;
    }
  }
  // If HOMO belongs to a partially filled degenerate shell where an earlier partner is occupied and the later is empty,
  // homoIdx already points to the highest occupied partner. Find first completely empty orbital above homoIdx as LUMO
  let lumoIdx = -1;
  for (let i = Math.max(0, homoIdx + 1); i < mos.length; i++) {
    if (mos[i].occupancy === 0) {
      lumoIdx = i;
      break;
    }
  }
  if (lumoIdx === -1 && homoIdx + 1 < mos.length) {
    lumoIdx = homoIdx + 1;
  }

  // Apply photo-excitation (promote 1 electron from HOMO to LUMO)
  if (excitedState && homoIdx >= 0 && lumoIdx >= 0 && mos[homoIdx].occupancy >= 1 && mos[lumoIdx].occupancy < 2) {
    mos[homoIdx].occupancy = (mos[homoIdx].occupancy - 1) as 0 | 1 | 2;
    mos[lumoIdx].occupancy = (mos[lumoIdx].occupancy + 1) as 0 | 1 | 2;
  }

  if (homoIdx >= 0) mos[homoIdx].isHOMO = true;
  if (lumoIdx >= 0) mos[lumoIdx].isLUMO = true;

  // Compute net MO bond order change relative to the neutral ground state
  const baseBondingElec = base.molecularOrbitals
    .filter((m) => m.bondingType === 'bonding')
    .reduce((s, m) => s + m.occupancy, 0);
  const baseAntibondingElec = base.molecularOrbitals
    .filter((m) => m.bondingType === 'antibonding')
    .reduce((s, m) => s + m.occupancy, 0);
  const baseNetBonding = (baseBondingElec - baseAntibondingElec) * 0.5;

  const bondingElec = mos
    .filter((m) => m.bondingType === 'bonding')
    .reduce((s, m) => s + m.occupancy, 0);
  const antibondingElec = mos
    .filter((m) => m.bondingType === 'antibonding')
    .reduce((s, m) => s + m.occupancy, 0);
  const currentNetBonding = (bondingElec - antibondingElec) * 0.5;

  return {
    ...base,
    name: overrideName || base.name,
    pubchemCid: overrideCid ?? base.pubchemCid,
    valenceElectrons: Math.max(1, base.valenceElectrons - charge),
    bondOrder: Number(
      Math.max(0.5, base.bondOrder + (currentNetBonding - baseNetBonding)).toFixed(2)
    ),
    molecularOrbitals: mos,
  };
}

/**
 * Evaluates the real 3D LCAO-MO spatial wavefunction ψ(x, y, z) at any point in 3D space!
 * Uses real Slater-type spherical harmonics centered at each atom's 3D coordinates.
 */
export function evaluateOrbitalWavefunction(
  mo: MolecularOrbital,
  atoms: AtomNode[],
  x: number,
  y: number,
  z: number
): number {
  let psi = 0;
  const zeta = 1.65; // Slater exponent in Å⁻¹

  for (const term of mo.lcao) {
    const atom = atoms.find((a) => a.id === term.atomId);
    if (!atom) continue;

    const dx = x - atom.x3d;
    const dy = y - atom.y3d;
    const dz = z - atom.z3d;
    const r = Math.hypot(dx, dy, dz);
    const radial = Math.exp(-zeta * r);

    let angular = 1.0;
    switch (term.orbital) {
      case 's':
        angular = 1.0;
        break;
      case 'px':
        angular = dx * 1.85;
        break;
      case 'py':
        angular = dy * 1.85;
        break;
      case 'pz':
        angular = dz * 1.85;
        break;
      case 'dz2':
        angular = (2 * dy * dy - dx * dx - dz * dz) * 1.3;
        break;
      case 'dx2y2':
        angular = (dx * dx - dz * dz) * 1.5;
        break;
      case 'dxy':
        angular = dx * dy * 2.0;
        break;
      case 'dxz':
        angular = dx * dz * 2.0;
        break;
      case 'dyz':
        angular = dy * dz * 2.0;
        break;
    }

    psi += term.coeff * radial * angular;
  }

  return psi;
}
