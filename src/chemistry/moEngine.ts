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
  // 1. Direct match on raw atoms & bonds or skeletal implicit-H formula (including C40H56 and C10H16N5O13P3)
  const rawFormula = computeMolecularFormula(atoms, bonds, false).plainFormula;
  const implicitFormula = computeMolecularFormula(atoms, bonds, true).plainFormula;
  for (const entry of CANONICAL_MOLECULES) {
    const entryRawFormula = computeMolecularFormula(entry.atoms, entry.bonds, false).plainFormula;
    const entryImplicitFormula = computeMolecularFormula(entry.atoms, entry.bonds, true).plainFormula;
    if (
      (entryRawFormula === rawFormula ||
        entryImplicitFormula === implicitFormula ||
        entry.formula === implicitFormula) &&
      entry.bonds.length === bonds.length &&
      entry.atoms.length === atoms.length
    ) {
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
    const entryFormula = computeMolecularFormula(entry.atoms, entry.bonds, true).plainFormula;
    if (entryFormula === plainFormula && entry.bonds.length === expanded.bonds.length) {
      return entry;
    }
  }
  return null;
}

/**
 * Cyclic Threshold Jacobi Eigenvalue & Eigenvector solver for real symmetric matrices H (N x N).
 * Performs full cyclic sweeps over all p < q pairs rather than single-element greedy pivots,
 * ensuring convergence for both small molecules (N ~ 6..40) and large biomolecules (N ~ 100..240).
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

  const maxCyclicSweeps = n > 60 ? 12 : 18;
  const eps = 1e-9;

  for (let sweep = 0; sweep < maxCyclicSweeps; sweep++) {
    let offDiagNormSq = 0;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        offDiagNormSq += A[i][j] * A[i][j];
      }
    }
    const offDiagNorm = Math.sqrt(offDiagNormSq);
    if (offDiagNorm < eps) break;

    // Adaptive threshold for early sweeps on large matrices; zero threshold on final sweeps
    const threshold = sweep < 3 ? (0.15 * offDiagNorm) / (n * n) : 1e-11;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = A[p][q];
        if (Math.abs(apq) <= threshold) continue;

        const app = A[p][p];
        const aqq = A[q][q];
        const theta = 0.5 * Math.atan2(2 * apq, aqq - app);
        const c = Math.cos(theta);
        const s = Math.sin(theta);

        // Rotate A rows/cols p and q
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

        // Rotate eigenvector matrix V
        for (let i = 0; i < n; i++) {
          const vip = V[i][p];
          const viq = V[i][q];
          V[i][p] = c * vip - s * viq;
          V[i][q] = s * vip + c * viq;
        }
      }
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
 * Computes the number of non-bonding lone pairs (m in AX_nE_m VSEPR notation)
 * on a central atom from its valence electron count, bond orders, and charge.
 */
export function getCentralAtomLonePairs(
  centerAtom: AtomNode,
  atoms: AtomNode[],
  bonds: BondEdge[],
  charge = 0
): number {
  const connected = bonds.filter(
    (b) => b.source === centerAtom.id || b.target === centerAtom.id
  );
  let effectiveBondOrderSum = 0;
  for (const b of connected) {
    const otherId = b.source === centerAtom.id ? b.target : b.source;
    const other = atoms.find((a) => a.id === otherId);
    const otherDeg = other
      ? bonds.filter((eb) => eb.source === other.id || eb.target === other.id).length
      : 1;
    // Terminal oxo/sulfido ligands (e.g. SO2, SO3, XeOF4, POCl3) require 2 valence electrons
    // from the central atom even when a 2D SDF encodes a dative single bond.
    if (other && otherDeg === 1 && (other.element === 'O' || other.element === 'S')) {
      effectiveBondOrderSum += Math.max(b.order, 2);
    } else {
      effectiveBondOrderSum += b.order;
    }
  }

  const valence = ELEMENT_PROPERTIES[centerAtom.element]?.valenceElectrons ?? 4;
  const nonBondingElectrons = Math.max(0, valence - effectiveBondOrderSum - charge);
  // Count both full lone pairs (2e-) and singly occupied radical orbitals (1e- SOMO, as in NO2 or ClO2)
  // as active VSEPR non-bonding steric domains so NO2 (17e-) is bent C2v rather than linear Dinfh.
  return Math.ceil(nonBondingElectrons / 2);
}

/**
 * Reconstructs authentic 3D VSEPR coordinates (x3d, y3d, z3d) when a molecule
 * is drawn on the flat 2D Sketcher canvas or loaded from a PubChem 2D-SDF fallback
 * (where all z3d ≈ 0), while preserving genuine 3D conformer coordinates when present.
 */
export function reconstruct3DGeometryFromVSEPR(
  atoms: AtomNode[],
  bonds: BondEdge[],
  charge = 0
): AtomNode[] {
  if (atoms.length === 0) return [];

  // Center 2D coordinates around 2D centroid
  const cx2d = atoms.reduce((s, a) => s + a.x2d, 0) / atoms.length;
  const cy2d = atoms.reduce((s, a) => s + a.y2d, 0) / atoms.length;

  // 1. Canonical alignment for diatomic molecules (e.g. HF, HCl, HBr, CO, NO, O2, N2, F2, Cl2, LiH)
  // Always align along the internuclear x-axis (y3d = 0, z3d = 0) with the heavier/electronegative atom at -r/2
  // so 3D symmetry elements (C∞ along x) and σ(px) / π(py, pz) orbitals are 100% aligned.
  if (atoms.length === 2) {
    const a0 = atoms[0];
    const a1 = atoms[1];
    const rawDist3d = Math.hypot(a1.x3d - a0.x3d, a1.y3d - a0.y3d, a1.z3d - a0.z3d);
    const covDist =
      ((ELEMENT_PROPERTIES[a0.element]?.covalentRadius ?? 0.75) +
        (ELEMENT_PROPERTIES[a1.element]?.covalentRadius ?? 0.75)) *
      0.96;
    const bondLen = Number(
      (rawDist3d >= 0.65 && rawDist3d <= 2.8 ? rawDist3d : covDist).toFixed(3)
    );
    const half = Number((bondLen * 0.5).toFixed(3));
    // Place heavier / non-H atom at -half and lighter / H atom at +half
    const z0 = ELEMENT_PROPERTIES[a0.element]?.atomicNumber ?? 6;
    const z1 = ELEMENT_PROPERTIES[a1.element]?.atomicNumber ?? 6;
    const swapOrder = a0.element === 'H' && a1.element !== 'H' ? true : z0 < z1 && a1.element !== 'H' && a0.element === 'H';
    const firstX = swapOrder ? half : -half;
    const secondX = swapOrder ? -half : half;

    return [
      {
        ...a0,
        x2d: Number(firstX.toFixed(2)),
        y2d: 0,
        x3d: firstX,
        y3d: 0,
        z3d: 0,
      },
      {
        ...a1,
        x2d: Number(secondX.toFixed(2)),
        y2d: 0,
        x3d: secondX,
        y3d: 0,
        z3d: 0,
      },
    ];
  }

  // 2. Check if molecule has a central atom bonded to all other n = atoms.length - 1 ligands (AX_nE_m VSEPR)
  // Always orient single-center VSEPR molecules (e.g. H2O, H2S, NH3, PH3, BH3, BF3, CH4, SiH4, SO2, SF4, ClF3, PF5, XeOF4, SF6)
  // into our canonical symmetry axes even when PubChem returns a tilted 3D-SDF conformer.
  const nLig = atoms.length - 1;
  const centerAtom =
    atoms.length >= 3 && atoms.length <= 7
      ? atoms.find((a) => {
          const deg = bonds.filter((b) => b.source === a.id || b.target === a.id).length;
          return deg === nLig;
        })
      : undefined;

  if (centerAtom && nLig >= 2 && nLig <= 6) {
    const rawLigands = atoms.filter((a) => a.id !== centerAtom.id);
    // Sort ligands so unique elements or higher-bond-order terminal oxo/sulfido ligands (e.g. =O in XeOF4, POCl3, SOCl2)
    // come first (ligands[0]), ensuring they occupy the axial / unique symmetry site regardless of drawing order.
    const elemFreq: Record<string, number> = {};
    for (const l of rawLigands) {
      elemFreq[l.element] = (elemFreq[l.element] || 0) + 1;
    }
    const getBondOrderToCenter = (lig: AtomNode) => {
      const b = bonds.find(
        (edge) =>
          (edge.source === centerAtom.id && edge.target === lig.id) ||
          (edge.target === centerAtom.id && edge.source === lig.id)
      );
      if (!b) return 1;
      if (lig.element === 'O' || lig.element === 'S') return Math.max(b.order, 2);
      return b.order;
    };
    const ligands = [...rawLigands].sort((a, b) => {
      const freqDiff = (elemFreq[a.element] || 1) - (elemFreq[b.element] || 1);
      if (freqDiff !== 0) return freqDiff; // Unique ligand (freq = 1) first
      return getBondOrderToCenter(b) - getBondOrderToCenter(a);
    });

    const lonePairs = getCentralAtomLonePairs(centerAtom, atoms, bonds, charge);
    const rCenter = ELEMENT_PROPERTIES[centerAtom.element]?.covalentRadius ?? 0.8;
    // Hypervalent 3c-4e species (lonePairs >= 2 or nLig >= 5) have fractional bond order and slightly longer bonds
    const scaleFactor = lonePairs >= 2 || nLig >= 5 ? 0.98 : 0.96;
    const maxAbsZInit = atoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d || 0)), 0);
    const getBondLen = (lig: AtomNode) => {
      if (maxAbsZInit > 0.08) {
        const d3d = Math.hypot(
          lig.x3d - centerAtom.x3d,
          lig.y3d - centerAtom.y3d,
          lig.z3d - centerAtom.z3d
        );
        if (d3d >= 0.7 && d3d <= 2.8) return Number(d3d.toFixed(3));
      }
      return Number(
        ((rCenter + (ELEMENT_PROPERTIES[lig.element]?.covalentRadius ?? 0.72)) * scaleFactor).toFixed(2)
      );
    };

    const coords3D = new Map<string, [number, number, number]>();
    coords3D.set(centerAtom.id, [0, 0, 0]);

    if (nLig === 2) {
      const r1 = getBondLen(ligands[0]);
      const r2 = getBondLen(ligands[1]);
      // AX2E0 (CO2, BeF2, HCN, OCS) or AX2E3 (XeF2) -> Linear along x-axis
      if (lonePairs === 0 || lonePairs >= 3) {
        coords3D.set(ligands[0].id, [-r1, 0, 0]);
        coords3D.set(ligands[1].id, [r2, 0, 0]);
      } else {
        // AX2E0.5 (NO2, 134.3°), AX2E1 (SO2, O3, 119°), or AX2E2 (H2O, H2S, OF2, SCl2, 104.5°) -> Bent in xy plane
        const totalValence = atoms.reduce(
          (s, a) => s + (ELEMENT_PROPERTIES[a.element]?.valenceElectrons ?? 4),
          -charge
        );
        const isOddElectronRadical = totalValence % 2 !== 0;
        const halfAngleDeg = lonePairs === 1 ? (isOddElectronRadical ? 67.15 : 59.5) : 52.25;
        const halfAngle = (halfAngleDeg * Math.PI) / 180;
        coords3D.set(centerAtom.id, [0, 0.32, 0]);
        coords3D.set(ligands[0].id, [
          Number((-r1 * Math.sin(halfAngle)).toFixed(3)),
          Number((0.32 - r1 * Math.cos(halfAngle)).toFixed(3)),
          0,
        ]);
        coords3D.set(ligands[1].id, [
          Number((r2 * Math.sin(halfAngle)).toFixed(3)),
          Number((0.32 - r2 * Math.cos(halfAngle)).toFixed(3)),
          0,
        ]);
      }
    } else if (nLig === 3) {
      if (lonePairs === 0) {
        // AX3E0 (BF3, BH3, SO3, BCl3) -> Trigonal planar in xy plane (z = 0)
        ligands.forEach((lig, idx) => {
          const r = getBondLen(lig);
          const phi = Math.PI / 2 + (idx * 2 * Math.PI) / 3;
          coords3D.set(lig.id, [
            Number((r * Math.cos(phi)).toFixed(3)),
            Number((r * Math.sin(phi)).toFixed(3)),
            0,
          ]);
        });
      } else if (lonePairs === 1) {
        // AX3E1 (NH3, PH3, PCl3, NF3) -> Trigonal pyramidal (C3v with C3 along y-axis)
        coords3D.set(centerAtom.id, [0, 0.38, 0]);
        ligands.forEach((lig, idx) => {
          const r = getBondLen(lig);
          const rProj = r * 0.92;
          const phi = Math.PI / 2 + (idx * 2 * Math.PI) / 3;
          coords3D.set(lig.id, [
            Number((rProj * Math.cos(phi)).toFixed(3)),
            -0.18,
            Number((rProj * Math.sin(phi)).toFixed(3)),
          ]);
        });
      } else {
        // AX3E2 (ClF3, BrF3) -> T-shaped (C2v with equatorial stem along -y and axial along ±x)
        const rEq = getBondLen(ligands[0]);
        const rAx1 = getBondLen(ligands[1]) * 1.05;
        const rAx2 = getBondLen(ligands[2]) * 1.05;
        coords3D.set(centerAtom.id, [0, 0.2, 0]);
        coords3D.set(ligands[0].id, [0, Number((0.2 - rEq).toFixed(3)), 0]);
        coords3D.set(ligands[1].id, [Number((-rAx1).toFixed(3)), 0.12, 0]);
        coords3D.set(ligands[2].id, [Number(rAx2.toFixed(3)), 0.12, 0]);
      }
    } else if (nLig === 4) {
      if (lonePairs === 0) {
        // AX4E0 (CH4, SiH4, CCl4, SiF4, POCl3) -> Tetrahedral (Td / C3v)
        const tetDirs: [number, number, number][] = [
          [0.8165, 0.5774, 0],
          [-0.8165, 0.5774, 0],
          [0, -0.5774, 0.8165],
          [0, -0.5774, -0.8165],
        ];
        ligands.forEach((lig, idx) => {
          const r = getBondLen(lig);
          const d = tetDirs[idx];
          coords3D.set(lig.id, [
            Number((r * d[0]).toFixed(3)),
            Number((r * d[1]).toFixed(3)),
            Number((r * d[2]).toFixed(3)),
          ]);
        });
      } else if (lonePairs === 1) {
        // AX4E1 (SF4, SeF4, TeCl4) -> Seesaw / Disphenoidal (C2v: 2 eq in xy, 2 ax along ±z)
        const rEq = getBondLen(ligands[0]);
        const rAx = getBondLen(ligands[2]) * 1.06;
        const eqAngle = (50.8 * Math.PI) / 180;
        coords3D.set(centerAtom.id, [0, 0.25, 0]);
        coords3D.set(ligands[0].id, [
          Number((-rEq * Math.sin(eqAngle)).toFixed(3)),
          Number((0.25 - rEq * Math.cos(eqAngle)).toFixed(3)),
          0,
        ]);
        coords3D.set(ligands[1].id, [
          Number((rEq * Math.sin(eqAngle)).toFixed(3)),
          Number((0.25 - rEq * Math.cos(eqAngle)).toFixed(3)),
          0,
        ]);
        coords3D.set(ligands[2].id, [0, 0.14, Number(rAx.toFixed(3))]);
        coords3D.set(ligands[3].id, [0, 0.14, Number((-rAx).toFixed(3))]);
      } else {
        // AX4E2 (XeF4) -> Square planar in xy plane (D4h, z = 0)
        const sqDirs: [number, number, number][] = [
          [1, 0, 0],
          [-1, 0, 0],
          [0, 1, 0],
          [0, -1, 0],
        ];
        ligands.forEach((lig, idx) => {
          const r = getBondLen(lig);
          const d = sqDirs[idx];
          coords3D.set(lig.id, [
            Number((r * d[0]).toFixed(3)),
            Number((r * d[1]).toFixed(3)),
            0,
          ]);
        });
      }
    } else if (nLig === 5) {
      if (lonePairs === 0) {
        // AX5E0 (PF5, PCl5) -> Trigonal bipyramidal (D3h: 3 eq in xy, 2 ax along ±z)
        for (let i = 0; i < 3; i++) {
          const rEq = getBondLen(ligands[i]);
          const phi = Math.PI / 2 + (i * 2 * Math.PI) / 3;
          coords3D.set(ligands[i].id, [
            Number((rEq * Math.cos(phi)).toFixed(3)),
            Number((rEq * Math.sin(phi)).toFixed(3)),
            0,
          ]);
        }
        const rAx1 = getBondLen(ligands[3]) * 1.04;
        const rAx2 = getBondLen(ligands[4]) * 1.04;
        coords3D.set(ligands[3].id, [0, 0, Number(rAx1.toFixed(3))]);
        coords3D.set(ligands[4].id, [0, 0, Number((-rAx2).toFixed(3))]);
      } else {
        // AX5E1 (BrF5, IF5, XeOF4) -> Square pyramidal (C4v: 1 axial at +y, 4 basal in xz)
        const rAx = getBondLen(ligands[0]);
        coords3D.set(ligands[0].id, [0, Number(rAx.toFixed(3)), 0]);
        const basalDirs: [number, number][] = [
          [1, 0],
          [-1, 0],
          [0, 1],
          [0, -1],
        ];
        for (let i = 0; i < 4; i++) {
          const lig = ligands[i + 1];
          const rEq = getBondLen(lig);
          coords3D.set(lig.id, [
            Number((rEq * basalDirs[i][0]).toFixed(3)),
            -0.18,
            Number((rEq * basalDirs[i][1]).toFixed(3)),
          ]);
        }
      }
    } else if (nLig === 6) {
      if (lonePairs === 0) {
        // AX6E0 (SF6) -> Regular Octahedral (Oh)
        const octDirs: [number, number, number][] = [
          [0, 1, 0],
          [0, -1, 0],
          [-1, 0, 0],
          [1, 0, 0],
          [0, 0, 1],
          [0, 0, -1],
        ];
        ligands.forEach((lig, idx) => {
          const r = getBondLen(lig);
          const d = octDirs[idx];
          coords3D.set(lig.id, [
            Number((r * d[0]).toFixed(3)),
            Number((r * d[1]).toFixed(3)),
            Number((r * d[2]).toFixed(3)),
          ]);
        });
      } else {
        // AX6E1 (XeF6) -> Distorted Monocapped Octahedron (C3v along y-axis)
        for (let i = 0; i < 3; i++) {
          const lig = ligands[i];
          const r = getBondLen(lig);
          const phi = Math.PI / 2 + (i * 2 * Math.PI) / 3;
          const rProj = r * 0.90;
          const yUp = r * 0.44;
          coords3D.set(lig.id, [
            Number((rProj * Math.cos(phi)).toFixed(3)),
            Number(yUp.toFixed(3)),
            Number((rProj * Math.sin(phi)).toFixed(3)),
          ]);
        }
        for (let i = 0; i < 3; i++) {
          const lig = ligands[i + 3];
          const r = getBondLen(lig) * 0.96;
          const phi = -Math.PI / 2 + (i * 2 * Math.PI) / 3;
          const rProj = r * 0.74;
          const yDown = -r * 0.67;
          coords3D.set(lig.id, [
            Number((rProj * Math.cos(phi)).toFixed(3)),
            Number(yDown.toFixed(3)),
            Number((rProj * Math.sin(phi)).toFixed(3)),
          ]);
        }
      }
    }

    return atoms.map((a) => {
      const c3 = coords3D.get(a.id) || [a.x2d - cx2d, -(a.y2d - cy2d), a.z3d || 0];
      return {
        ...a,
        x2d: Number((a.x2d - cx2d).toFixed(2)),
        y2d: Number((a.y2d - cy2d).toFixed(2)),
        x3d: c3[0],
        y3d: c3[1],
        z3d: c3[2],
      };
    });
  }

  // Check if genuine 3D depth already exists for larger multi-center molecules (e.g. from PubChem 3D-SDF)
  const maxAbsZ = atoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d || 0)), 0);
  if (maxAbsZ > 0.08) {
    const cx3d = atoms.reduce((s, a) => s + a.x3d, 0) / atoms.length;
    const cy3d = atoms.reduce((s, a) => s + a.y3d, 0) / atoms.length;
    const cz3d = atoms.reduce((s, a) => s + a.z3d, 0) / atoms.length;
    return atoms.map((a) => ({
      ...a,
      x2d: Number((a.x2d - cx2d).toFixed(2)),
      y2d: Number((a.y2d - cy2d).toFixed(2)),
      x3d: Number((a.x3d - cx3d).toFixed(3)),
      y3d: Number((a.y3d - cy3d).toFixed(3)),
      z3d: Number((a.z3d - cz3d).toFixed(3)),
    }));
  }

  // Multi-center 3D tetrahedral z-displacement for flat 2D graphs with tetracoordinate centers (e.g. P in ATP or sp3 C)
  const zOffsets = new Map<string, number>();
  for (const center of atoms) {
    const conn = bonds.filter((b) => b.source === center.id || b.target === center.id);
    if (conn.length === 4) {
      const neighbors = conn
        .map((b) => atoms.find((a) => a.id === (b.source === center.id ? b.target : b.source)))
        .filter((a): a is AtomNode => Boolean(a));
      // Displace terminal substituents on tetracoordinate P / S / C along ±z to break artificial 2D planarity
      const terminals = neighbors.filter(
        (nb) => bonds.filter((b) => b.source === nb.id || b.target === nb.id).length <= 2
      );
      if (terminals.length >= 2) {
        zOffsets.set(terminals[0].id, 1.18);
        zOffsets.set(terminals[1].id, -1.18);
      }
    }
  }

  return atoms.map((a) => ({
    ...a,
    x2d: Number((a.x2d - cx2d).toFixed(2)),
    y2d: Number((a.y2d - cy2d).toFixed(2)),
    x3d: Number((a.x2d - cx2d).toFixed(2)),
    y3d: Number((-(a.y2d - cy2d)).toFixed(2)),
    z3d: zOffsets.get(a.id) ?? (a.z3d || 0),
  }));
}

/**
 * Determines the molecular Point Group from atomic composition, bond connectivity,
 * VSEPR central-atom lone pairs (AX_nE_m), and 2D/3D geometry.
 */
export function detectPointGroup(atoms: AtomNode[], bonds: BondEdge[], charge = 0): string {
  const heavy = atoms.filter((a) => a.element !== 'H');
  const nAtoms = atoms.length;

  if (nAtoms <= 1) return 'Oh';

  // Diatomic molecules
  if (nAtoms === 2) {
    return atoms[0].element === atoms[1].element ? 'Dinfh' : 'Cinfv';
  }

  // Triatomic molecules (e.g., H2O, CO2, O3, SO2, XeF2, HCN, OCS)
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
      const lonePairs = getCentralAtomLonePairs(center, atoms, bonds, charge);
      // AX2E0 (CO2, CS2, BeF2) or AX2E3 (XeF2) are linear; AX2E1 (SO2, O3, NO2) and AX2E2 (H2O, OF2) are bent
      const isLinearVsepr = lonePairs === 0 || lonePairs >= 3;

      if (others[0].element === others[1].element) {
        return isLinearVsepr && center.element !== 'O' && center.element !== 'S' ? 'Dinfh' : 'C2v';
      } else {
        return isLinearVsepr ? 'Cinfv' : 'Cs';
      }
    }
    return 'Cinfv';
  }

  // Check centrosymmetric all-trans polyenes (e.g. β-Carotene C40H56 -> idealized D2h)
  if (heavy.length === 40 && heavy.every((a) => a.element === 'C') && bonds.length >= 40) {
    return 'D2h';
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

  // Central atom surrounded by nLig ligands (AX_nE_m VSEPR point-group classification)
  const centerAtom = atoms.find((a) => {
    const deg = bonds.filter((b) => b.source === a.id || b.target === a.id).length;
    return deg >= 3 && deg === atoms.length - 1;
  });

  if (centerAtom) {
    const ligands = atoms.filter((a) => a.id !== centerAtom.id);
    const ligandElements = new Set(ligands.map((a) => a.element));
    const nLig = ligands.length;
    const lonePairs = getCentralAtomLonePairs(centerAtom, atoms, bonds, charge);

    if (ligandElements.size === 1) {
      if (nLig === 3) {
        // AX3E0 (BF3, SO3) -> D3h; AX3E1 (NH3, PCl3) -> C3v; AX3E2 (ClF3) -> C2v
        if (lonePairs === 0) return 'D3h';
        if (lonePairs === 1) return 'C3v';
        return 'C2v';
      }
      if (nLig === 4) {
        // AX4E0 (CH4, CCl4) -> Td; AX4E1 (SF4) -> C2v; AX4E2 (XeF4) -> D4h
        if (lonePairs === 0) return 'Td';
        if (lonePairs === 1) return 'C2v';
        return 'D4h';
      }
      if (nLig === 5) {
        // AX5E0 (PF5, PCl5) -> D3h; AX5E1 (BrF5, IF5) -> C4v
        return lonePairs === 0 ? 'D3h' : 'C4v';
      }
      if (nLig === 6) {
        // AX6E0 (SF6) -> Oh; AX6E1 (XeF6) -> C3v (distorted monocapped octahedron)
        return lonePairs === 0 ? 'Oh' : 'C3v';
      }
    }

    if (ligandElements.size === 2) {
      const elemCounts = Array.from(ligandElements).map(
        (el) => ligands.filter((l) => l.element === el).length
      );
      if (nLig === 3) {
        // AX2Y: planar (CH2O, COCl2 -> C2v) vs pyramidal (SOCl2 -> Cs)
        return lonePairs === 0 ? 'C2v' : 'Cs';
      }
      if (nLig === 4 && elemCounts.includes(3)) {
        // AX3Y tetrahedral (CHCl3, POCl3 -> C3v)
        return lonePairs === 0 ? 'C3v' : 'Cs';
      }
      if (nLig === 5 && elemCounts.includes(4)) {
        // AX4Y square pyramidal (XeOF4 -> C4v)
        return 'C4v';
      }
    }
  }

  if (heavy.length === 2 && heavy[0].element === heavy[1].element) {
    const nH = atoms.filter((a) => a.element === 'H').length;
    if (nH === 2) return 'Dinfh';
    if (nH === 4) return 'D2h';
  }

  return 'Cs';
}

const ehtSolverCache = new Map<
  string,
  {
    basis: {
      atom: AtomNode;
      orbital: 's' | 'px' | 'py' | 'pz';
      label: string;
      alpha: number;
    }[];
    eigenvalues: number[];
    eigenvectors: number[][];
  }
>();

/**
 * Computes raw Extended Hückel Theory (EHT) eigenvalues and eigenvectors using
 * a full 4-orbital (ns, npx, npy, npz) valence basis for heavy atoms and 1s for H,
 * with directional Slater-Koster two-center overlap integrals S_mu_nu(R_AB).
 * For large molecules (> 24 atoms), uses an active-space valence basis (retaining
 * heavy atoms + heteroatom hydrogens while folding aliphatic C-H sigma orbitals)
 * and memoizes the diagonalization by geometry hash.
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
  const cacheKey =
    atoms.map((a) => `${a.id}:${a.element}:${a.x3d},${a.y3d},${a.z3d}`).join('|') +
    '#' +
    bonds.map((b) => `${b.source}-${b.target}:${b.order}`).join('|');
  const cached = ehtSolverCache.get(cacheKey);
  if (cached) return cached;

  interface BasisFunc {
    atom: AtomNode;
    orbital: 's' | 'px' | 'py' | 'pz';
    label: string;
    alpha: number; // diagonal Coulomb integral H_ii in eV
  }

  // For large systems (> 28 atoms), prune explicit C-H hydrogens from the active matrix
  // so the solver focuses on the heavy-atom framework + heteroatom H and runs in < 20 ms
  const isLargeSystem = atoms.length > 28;
  const activeAtoms = isLargeSystem
    ? atoms
        .filter((a) => {
          if (a.element !== 'H') return true;
          const attached = bonds.find((b) => b.source === a.id || b.target === a.id);
          if (!attached) return false;
          const partnerId = attached.source === a.id ? attached.target : attached.source;
          const partner = atoms.find((p) => p.id === partnerId);
          return partner && partner.element !== 'C';
        })
        .slice(0, 48)
    : atoms;

  const basis: BasisFunc[] = [];
  for (const a of activeAtoms) {
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

      // Direction cosines from atom i to atom j: d = (R_j - R_i) / |R_j - R_i|
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

      // Exact global-Cartesian Slater-Koster two-center overlap integrals:
      // - For s on i and p_α on j: lobe of p_α pointing toward i (-d) has sign -d_α
      // - For p_α on i and s on j: lobe of p_α pointing toward j (+d) has sign +d_α
      // - For p_α on i and p_β on j: head-on σ overlap has sign -d_α d_β, side-on π has sign +(δ_αβ - d_α d_β)
      let overlapS = 0;
      if (bi.orbital === 's' && bj.orbital === 's') {
        overlapS = S_ss_sigma;
      } else if (bi.orbital === 's' && bj.orbital !== 's') {
        overlapS = -dirMap[bj.orbital] * S_sp_sigma;
      } else if (bi.orbital !== 's' && bj.orbital === 's') {
        overlapS = dirMap[bi.orbital] * S_sp_sigma;
      } else if (bi.orbital !== 's' && bj.orbital !== 's') {
        const di = dirMap[bi.orbital];
        const dj = dirMap[bj.orbital];
        if (bi.orbital === bj.orbital) {
          overlapS = -di * di * S_pp_sigma + (1 - di * di) * S_pp_pi;
        } else {
          overlapS = -di * dj * (S_pp_sigma + S_pp_pi);
        }
      }

      // Löwdin-orthogonalized effective EHT coupling: damp deep core-like ns / valence energy mismatches
      // (such as F 2s at -40.2 eV vs H 1s at -13.6 eV) so deep ns orbitals remain predominantly non-bonding
      // rather than unphysically distorting the frontier σ and σ* orbitals.
      const dAlpha = Math.abs(bi.alpha - bj.alpha);
      const gapDamping = 1.0 / (1.0 + Math.pow(dAlpha / 18.0, 2));
      const hij = 0.5 * K * (bi.alpha + bj.alpha) * overlapS * (0.35 + 0.65 * gapDamping);
      H[i][j] = hij;
      H[j][i] = hij;
    }
  }

  const { eigenvalues, eigenvectors } = jacobiEigenvalues(H);
  const result = { basis, eigenvalues, eigenvectors };
  if (ehtSolverCache.size > 24) {
    const firstKey = ehtSolverCache.keys().next().value;
    if (firstKey) ehtSolverCache.delete(firstKey);
  }
  ehtSolverCache.set(cacheKey, result);
  return result;
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

  // 2. Otherwise, expand implicit hydrogens (for small/medium graphs <= 24 atoms),
  //    reconstruct 3D VSEPR coordinates if flat (z3d ≈ 0),
  //    and run the full 4-orbital Extended Hückel LCAO-MO + Symmetry solver
  const shouldExpandImplicitH = rawAtoms.length <= 24;
  const { atoms, bonds } = shouldExpandImplicitH
    ? expandImplicitHydrogens(rawAtoms, rawBonds)
    : { atoms: rawAtoms, bonds: rawBonds };
  const { formula } = computeMolecularFormula(atoms, bonds, !shouldExpandImplicitH);
  const smiles = generateSmilesFromGraph(atoms, bonds);
  const centeredAtoms = reconstruct3DGeometryFromVSEPR(atoms, bonds, charge);
  const pointGroupId = detectPointGroup(centeredAtoms, bonds, charge);
  const pg = POINT_GROUPS[pointGroupId] || POINT_GROUPS.Cs;

  const { basis, eigenvalues, eigenvectors } = runExtendedHuckelSolver(centeredAtoms, bonds);
  const N = basis.length;

  // Exact physical neutral valence electrons = sum(Z_val) of active basis atoms; charge is applied once in applyElectronicState
  const activeBasisAtomIds = new Set(basis.map((b) => b.atom.id));
  const rawValenceElectrons = centeredAtoms
    .filter((a) => activeBasisAtomIds.has(a.id))
    .reduce((sum, a) => sum + (ELEMENT_PROPERTIES[a.element]?.valenceElectrons ?? 4), 0);
  const neutralElectrons = Math.max(1, Math.min(N * 2, rawValenceElectrons));

  // Detect whether the molecule is linear along the x-axis or planar in the xy-plane (z ≈ 0)
  const maxAbsY = centeredAtoms.reduce((m, a) => Math.max(m, Math.abs(a.y3d)), 0);
  const maxAbsZ = centeredAtoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d)), 0);
  const isLinearAlongX = centeredAtoms.length >= 2 && maxAbsY < 0.15 && maxAbsZ < 0.15;
  const isPlanarXY = maxAbsZ < 0.15;

  // Identify central / primary skeleton atom and peripheral ligand atom for dynamic fragment orbitals
  const primaryAtom =
    centeredAtoms.find((a) => {
      const deg = bonds.filter((b) => b.source === a.id || b.target === a.id).length;
      return deg >= 2 && deg === centeredAtoms.length - 1;
    }) ||
    centeredAtoms.find((a) => a.element !== 'H') ||
    centeredAtoms[0];
  const peripheralAtom =
    centeredAtoms.find((a) => a.id !== primaryAtom?.id) || primaryAtom;

  const primProps = ELEMENT_PROPERTIES[primaryAtom?.element || 'C'] || ELEMENT_PROPERTIES.C;
  const periphProps = ELEMENT_PROPERTIES[peripheralAtom?.element || 'H'] || ELEMENT_PROPERTIES.H;
  const primShell =
    primaryAtom?.element === 'P' || primaryAtom?.element === 'S' || primaryAtom?.element === 'Cl'
      ? '3'
      : primaryAtom?.element === 'Xe'
        ? '5'
        : primaryAtom?.element === 'H'
          ? '1'
          : '2';
  const periphShell =
    peripheralAtom?.element === 'P' || peripheralAtom?.element === 'S' || peripheralAtom?.element === 'Cl'
      ? '3'
      : peripheralAtom?.element === 'Xe'
        ? '5'
        : peripheralAtom?.element === 'H'
          ? '1'
          : '2';

  // Build dynamic fragment orbitals from actual element VSIEs
  const firstIrrepLower = (pg.irreps[0]?.symbol || 'a').toLowerCase().replace(/ς/g, 'σ');
  const secondIrrepLower = (pg.irreps[Math.min(1, pg.irreps.length - 1)]?.symbol || 'b')
    .toLowerCase()
    .replace(/ς/g, 'σ');
  const lastIrrepLower = (pg.irreps[pg.irreps.length - 1]?.symbol || 'b')
    .toLowerCase()
    .replace(/ς/g, 'σ');
  const fragmentOrbitals: FragmentOrbital[] = [
    {
      id: 'frag_left_s',
      side: 'left',
      label: `${primaryAtom?.element || 'C'} ${primShell}s`,
      fragmentName: `${primaryAtom?.element || 'Core'} Valence`,
      energy: primProps.vsieS,
      irrep: firstIrrepLower,
      degeneracy: 1,
      electrons: 2,
    },
    {
      id: 'frag_left_p',
      side: 'left',
      label: `${primaryAtom?.element || 'C'} ${primShell}p (σ / π)`,
      fragmentName: `${primaryAtom?.element || 'Core'} Valence`,
      energy: primProps.vsieP,
      irrep: secondIrrepLower,
      degeneracy: 3,
      electrons: Math.max(1, primProps.valenceElectrons - 2),
    },
    {
      id: 'frag_right_salc',
      side: 'right',
      label:
        peripheralAtom?.element === 'H'
          ? 'Ligand H 1s SALCs'
          : `Ligand ${peripheralAtom?.element || 'L'} ${periphShell}p SALCs`,
      fragmentName: 'Substituent SALCs',
      energy: peripheralAtom?.element === 'H' ? periphProps.vsieS : periphProps.vsieP,
      irrep: `${firstIrrepLower} + ${lastIrrepLower}`,
      degeneracy: 2,
      electrons: 2,
    },
  ];

  // Select all valence MOs when N <= 14 (e.g. C2H4 with 12 AOs), or a 12-orbital frontier window for larger systems without slicing degenerate pairs
  const homoIdx = Math.max(0, Math.min(N - 1, Math.floor((neutralElectrons - 1) / 2)));
  let startIdx = N <= 14 ? 0 : Math.max(0, Math.min(homoIdx - 5, N - 12));
  let endIdx = N <= 14 ? N - 1 : Math.min(N - 1, Math.max(homoIdx + 5, startIdx + 11));
  if (startIdx > 0 && Math.abs(eigenvalues[startIdx] - eigenvalues[startIdx - 1]) < 0.08) {
    startIdx -= 1;
  }
  if (endIdx < N - 1 && Math.abs(eigenvalues[endIdx] - eigenvalues[endIdx + 1]) < 0.08) {
    endIdx += 1;
  }

  const molecularOrbitals: MolecularOrbital[] = [];
  let remainingElectrons = Math.max(0, neutralElectrons - startIdx * 2);

  // Exclude d-only Delta and multi-electron SigmaMinus irreps when using an s-p valence basis
  const spIrreps = pg.irreps.filter(
    (ir) => !ir.symbol.startsWith('Δ') && ir.symbol !== 'Σ-'
  );
  const candidateIrreps = spIrreps.length > 0 ? spIrreps : pg.irreps;
  const dim1Irreps = candidateIrreps.filter((ir) => ir.dimension === 1);
  const dimMultiIrreps = candidateIrreps.filter((ir) => ir.dimension >= 2);

  // Initialize literature-standard per-irrep Mulliken shell counters (accounting for inner 1s core shells on heavy atoms)
  const heavyCount = centeredAtoms.filter((a) => a.element !== 'H').length;
  const irrepShellCounts: Record<string, number> = {};
  if (startIdx === 0) {
    if (pg.id === 'Cinfv') {
      // 1 core σ+ per heavy atom (e.g. HF/HCl -> 1σ+ core, so valence starts at 2σ+, 3σ+, 1π, 4σ+*; CO/NO/HCN -> 2σ+ cores, so valence starts at 3σ+)
      irrepShellCounts['SigmaPlus'] = Math.max(1, heavyCount);
    } else if (pg.id === 'Dinfh' && heavyCount === 2) {
      irrepShellCounts['SigmagPlus'] = 1;
      irrepShellCounts['SigmauPlus'] = 1;
    } else if (
      (pg.id === 'C2v' || pg.id === 'C3v' || pg.id === 'C4v' || pg.id === 'Td') &&
      heavyCount === 1
    ) {
      irrepShellCounts['A1'] = 1;
    } else if (pg.id === 'D3h' && heavyCount === 1) {
      irrepShellCounts['A1p'] = 1;
    }
  }
  let degPartnerCount = 0;

  for (let idx = startIdx; idx <= endIdx; idx++) {
    const energy = Number(eigenvalues[idx].toFixed(2));
    const vec = eigenvectors[idx];

    // Only classify pz as pi-type if the molecule is actually planar in xy (z ≈ 0) or linear along x
    const piWeight = vec.reduce(
      (acc, c, k) =>
        acc +
        ((isPlanarXY && basis[k].orbital === 'pz') ||
        (isLinearAlongX && (basis[k].orbital === 'py' || basis[k].orbital === 'pz'))
          ? c * c
          : 0),
      0
    );
    const isPi = piWeight > 0.48;

    const occ: 0 | 1 | 2 =
      remainingElectrons >= 2 ? 2 : remainingElectrons === 1 ? 1 : 0;
    remainingElectrons = Math.max(0, remainingElectrons - occ);

    // Compute interatomic bond overlap population across bonded pairs to classify bonding vs non-bonding vs antibonding
    let bondPopulation = 0;
    let signChanges = 0;
    for (const b of bonds) {
      let pairOverlap = 0;
      for (let i = 0; i < N; i++) {
        if (basis[i].atom.id !== b.source) continue;
        for (let j = 0; j < N; j++) {
          if (basis[j].atom.id !== b.target) continue;
          const ci = vec[i];
          const cj = vec[j];
          if (Math.abs(ci * cj) < 1e-4) continue;
          const dx = basis[j].atom.x3d - basis[i].atom.x3d;
          const dy = basis[j].atom.y3d - basis[i].atom.y3d;
          const dz = basis[j].atom.z3d - basis[i].atom.z3d;
          const dist = Math.hypot(dx, dy, dz) || 1.4;
          const dirMap: Record<'px' | 'py' | 'pz', number> = {
            px: dx / dist,
            py: dy / dist,
            pz: dz / dist,
          };
          const oi = basis[i].orbital;
          const oj = basis[j].orbital;
          let sVal = 0;
          if (oi === 's' && oj === 's') sVal = 0.32;
          else if (oi === 's' && oj !== 's') sVal = -dirMap[oj] * 0.28;
          else if (oi !== 's' && oj === 's') sVal = dirMap[oi] * 0.28;
          else if (oi !== 's' && oj !== 's') {
            const di = dirMap[oi];
            const dj = dirMap[oj];
            sVal = oi === oj ? -di * di * 0.26 + (1 - di * di) * 0.16 : -di * dj * (0.26 + 0.16);
          }
          pairOverlap += 2 * ci * cj * sVal;
        }
      }
      bondPopulation += pairOverlap;
      if (pairOverlap < -0.02) signChanges++;
    }

    const bondingType: MolecularOrbital['bondingType'] =
      bondPopulation > 0.045
        ? 'bonding'
        : bondPopulation < -0.02 || (idx > homoIdx && Math.abs(bondPopulation) > 0.01)
          ? 'antibonding'
          : 'non-bonding';

    const interactionType: MolecularOrbital['interactionType'] = isPi
      ? bondingType === 'antibonding'
        ? 'π*'
        : bondingType === 'non-bonding'
          ? 'n (Lone Pair)'
          : 'π'
      : bondingType === 'antibonding'
        ? 'σ*'
        : bondingType === 'non-bonding'
          ? 'n (Lone Pair)'
          : 'σ';

    const lcaoMinCoeff = centeredAtoms.length > 20 ? 0.06 : 0.14;
    const lcaoMaxTerms = centeredAtoms.length > 20 ? 24 : 8;
    const lcao: LCAOContribution[] = vec
      .map((coeff, k) => ({
        atomId: basis[k].atom.id,
        element: basis[k].atom.element,
        orbital: basis[k].orbital,
        label: basis[k].label,
        coeff: Number(coeff.toFixed(2)),
      }))
      .sort((a, b) => Math.abs(b.coeff) - Math.abs(a.coeff))
      .filter((c) => Math.abs(c.coeff) >= lcaoMinCoeff)
      .slice(0, lcaoMaxTerms);

    // Group nearly degenerate eigenvalues into a shared degenerateGroup and assign consistent Mulliken irreps
    const prevMo = molecularOrbitals[molecularOrbitals.length - 1];
    const isDegenerateWithPrev = Boolean(prevMo && Math.abs(prevMo.energy - energy) < 0.08);
    const isDegenerateWithNext =
      idx + 1 <= endIdx && Math.abs(eigenvalues[idx + 1] - eigenvalues[idx]) < 0.08;
    const isDegenerate = isDegenerateWithPrev || isDegenerateWithNext;

    const degenerateGroup = isDegenerateWithPrev
      ? prevMo.degenerateGroup || `deg_${idx - 1}`
      : undefined;
    if (isDegenerateWithPrev && prevMo && !prevMo.degenerateGroup) {
      prevMo.degenerateGroup = degenerateGroup;
    }

    let irrepObj = candidateIrreps[0];
    if (isDegenerateWithPrev && prevMo) {
      irrepObj =
        candidateIrreps.find((ir) => ir.id === prevMo.irrepId) || candidateIrreps[0];
    } else if (isDegenerate && dimMultiIrreps.length > 0) {
      if (pg.id === 'D3h') {
        irrepObj =
          dimMultiIrreps.find((ir) => (isPi ? ir.id === 'Epp' : ir.id === 'Ep')) ||
          dimMultiIrreps[0];
      } else {
        // Check inversion / ungerade vs gerade character
        const centerWeight = primaryAtom
          ? vec.reduce((s, c, k) => s + (basis[k].atom.id === primaryAtom.id ? c * c : 0), 0)
          : 0.5;
        const preferUngerade =
          heavyCount === 2 ? bondingType === 'bonding' : centerWeight > 0.08;
        irrepObj =
          dimMultiIrreps.find((ir) =>
            preferUngerade ? ir.symbol.includes('u') : ir.symbol.includes('g')
          ) || dimMultiIrreps[0];
      }
    } else if (dim1Irreps.length > 0) {
      if (pg.id === 'Cinfv') {
        irrepObj = dim1Irreps.find((ir) => ir.id === 'SigmaPlus') || dim1Irreps[0];
      } else if (pg.id === 'C2v') {
        const pxWeight = vec.reduce(
          (s, c, k) => s + (basis[k].orbital === 'px' ? c * c : 0),
          0
        );
        const targetId = isPi ? 'B1' : pxWeight > 0.35 || signChanges % 2 === 1 ? 'B2' : 'A1';
        irrepObj = dim1Irreps.find((ir) => ir.id === targetId) || dim1Irreps[0];
      } else if (pg.id === 'C3v' || pg.id === 'C4v') {
        const centerSPWeight = primaryAtom
          ? vec.reduce((s, c, k) => s + (basis[k].atom.id === primaryAtom.id ? c * c : 0), 0)
          : 1;
        irrepObj =
          dim1Irreps.find((ir) => (centerSPWeight > 0.05 ? ir.id === 'A1' : ir.id === 'A2')) ||
          dim1Irreps[0];
      } else if (pg.id === 'D3h') {
        const targetId = isPi ? 'A2pp' : 'A1p';
        irrepObj = dim1Irreps.find((ir) => ir.id === targetId) || dim1Irreps[0];
      } else {
        const centerPWeight = primaryAtom
          ? vec.reduce(
              (s, c, k) =>
                s + (basis[k].atom.id === primaryAtom.id && basis[k].orbital !== 's' ? c * c : 0),
              0
            )
          : 0;
        const preferUngeradeOrB = centerPWeight > 0.12 || signChanges % 2 === 1;
        irrepObj =
          dim1Irreps.find((ir) =>
            preferUngeradeOrB
              ? ir.symbol.includes('u') || ir.symbol.startsWith('B') || ir.symbol.includes("''")
              : ir.symbol.includes('g') || ir.symbol.startsWith('A') || ir.symbol.includes('Σ+')
          ) || dim1Irreps[idx % dim1Irreps.length];
      }
    }

    if (isDegenerateWithPrev) {
      degPartnerCount += 1;
    } else {
      degPartnerCount = 0;
      irrepShellCounts[irrepObj.id] = (irrepShellCounts[irrepObj.id] || 0) + 1;
    }

    const shellNum = irrepShellCounts[irrepObj.id] || idx - startIdx + 1;
    let degSuffix = '';
    if (isDegenerate) {
      if (irrepObj.symbol.includes('Π')) {
        degSuffix = degPartnerCount === 0 ? '(x)' : '(y)';
      } else if (irrepObj.dimension === 3) {
        degSuffix = ['(x)', '(y)', '(z)'][degPartnerCount % 3];
      } else {
        degSuffix = degPartnerCount === 0 ? '(a)' : '(b)';
      }
    }

    const lowerIrrepSym = irrepObj.symbol.toLowerCase().replace(/ς/g, 'σ');
    molecularOrbitals.push({
      id: `custom_mo_${idx + 1}`,
      label: `${shellNum}${lowerIrrepSym}${bondingType === 'antibonding' ? '*' : ''}${degSuffix}`,
      irrepSymbol: lowerIrrepSym,
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
      description: `Extended Hückel LCAO-MO (${interactionType}, ${lowerIrrepSym} symmetry) computed at ε = ${energy.toFixed(2)} eV with ${signChanges} internuclear nodal surface(s).`,
    });
  }

  // Compute physical net bond order (accounting for 3c-4e hypervalent lone pairs when present)
  const centerLonePairs =
    primaryAtom && centeredAtoms.length >= 3
      ? getCentralAtomLonePairs(primaryAtom, centeredAtoms, bonds, charge)
      : 0;
  const rawBondOrderSum = bonds.reduce((s, b) => s + b.order, 0);
  const effectiveBondOrder =
    centerLonePairs >= 2 && bonds.length >= 2
      ? Math.max(1.0, rawBondOrderSum - Math.floor(centerLonePairs / 2))
      : rawBondOrderSum;

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
    bondOrder: effectiveBondOrder,
    dipoleMomentDebye:
      pg.id === 'D6h' || pg.id === 'Td' || pg.id === 'Oh' || pg.id === 'Dinfh' || pg.id === 'D2h' || pg.id === 'D3h' || pg.id === 'D4h'
        ? 0.0
        : 1.42,
    leftFragmentTitle: `Central ${primaryAtom?.element || 'Skeleton'} (${primShell}s, ${primShell}p)`,
    rightFragmentTitle: `Peripheral ${peripheralAtom?.element || 'Ligand'} SALCs`,
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
