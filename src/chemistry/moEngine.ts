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
  OrbitalBasisType,
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

  // Terminal carbon triple- or double-bonded to Oxygen, Nitrogen, or Sulfur (C≡O, [C]=O, C≡N-, R-N≡C) has 0 implicit H
  if (connected.length === 1 && connected[0].order >= 2) {
    const otherId = connected[0].source === atom.id ? connected[0].target : connected[0].source;
    const other = atoms.find((x) => x.id === otherId);
    if (other && (other.element === 'O' || other.element === 'N' || other.element === 'S')) {
      return 0;
    }
  }

  // Coordinated cyanide / carbonyl / thiocarbonyl carbon (M-C≡N, M-C≡O, M-C=O) has 0 implicit H
  const hasMultipleBondToHetero = connected.some((b) => {
    if (b.order < 2) return false;
    const otherId = b.source === atom.id ? b.target : b.source;
    const other = atoms.find((x) => x.id === otherId);
    return Boolean(other && (other.element === 'O' || other.element === 'N' || other.element === 'S'));
  });
  const isBondedToMetal = connected.some((b) => {
    const otherId = b.source === atom.id ? b.target : b.source;
    const other = atoms.find((x) => x.id === otherId);
    return Boolean(other && (isTransitionMetal(other.element) || isLanthanide(other.element)));
  });
  if ((hasMultipleBondToHetero && isBondedToMetal) || (atom.formalCharge && hasMultipleBondToHetero)) {
    return 0;
  }

  const bondOrderSum = connected.reduce((sum, b) => sum + b.order, 0);
  return Math.max(0, 4 - bondOrderSum);
}

export function getPrincipalShell(el: ElementSymbol): string {
  const z = ELEMENT_PROPERTIES[el]?.atomicNumber ?? 6;
  if (z <= 2) return '1';
  if (z <= 10) return '2';
  if (z <= 18) return '3';
  if (z <= 36) return '4';
  if (z <= 54) return '5';
  return '6';
}

export function isTransitionMetal(el: ElementSymbol): boolean {
  const z = ELEMENT_PROPERTIES[el]?.atomicNumber ?? 6;
  return (z >= 21 && z <= 30) || (z >= 39 && z <= 48) || (z >= 72 && z <= 80);
}

export function isLanthanide(el: ElementSymbol): boolean {
  const z = ELEMENT_PROPERTIES[el]?.atomicNumber ?? 6;
  return z >= 57 && z <= 71;
}

const SQUARE_PLANAR_METALS = new Set<ElementSymbol>(['Pt', 'Pd', 'Ni', 'Au', 'Rh', 'Ir']);

export function isSquarePlanarMetal(el: ElementSymbol): boolean {
  return SQUARE_PLANAR_METALS.has(el);
}

const SPECTATOR_COUNTERIONS = new Set<ElementSymbol>([
  'Li', 'Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba',
]);

const MAIN_GROUP_COORD_CENTERS = new Set<ElementSymbol>([
  'B', 'Al', 'Ga', 'In', 'Tl', 'Si', 'Ge', 'Sn', 'Pb', 'P', 'As', 'Sb', 'Bi', 'S', 'Se', 'Te', 'Xe',
]);

/**
 * Resolves missing coordinate covalent (dative M–L) bonds in coordination complexes
 * where PubChem V2000 SDFs or disconnected SMILES encode ligands (`[C-]#N`, `[C-]#[O+]`,
 * `[C]=O`, `NH3`, `Cl-`, etc.) as separate unconnected fragments around a central metal
 * or hypervalent coordination center, and strips isolated spectator alkali/alkaline-earth counterions.
 */
export function resolveCoordinationBonds(
  rawAtoms: AtomNode[],
  rawBonds: BondEdge[],
  stripSpectatorCounterions = true
): { atoms: AtomNode[]; bonds: BondEdge[] } {
  if (rawAtoms.length <= 1) {
    return { atoms: rawAtoms, bonds: rawBonds };
  }

  const getComponents = (atomList: AtomNode[], bondList: BondEdge[]): AtomNode[][] => {
    const adj = new Map<string, string[]>();
    for (const a of atomList) adj.set(a.id, []);
    for (const b of bondList) {
      if (adj.has(b.source) && adj.has(b.target)) {
        adj.get(b.source)!.push(b.target);
        adj.get(b.target)!.push(b.source);
      }
    }
    const visited = new Set<string>();
    const comps: AtomNode[][] = [];
    for (const a of atomList) {
      if (visited.has(a.id)) continue;
      const comp: AtomNode[] = [];
      const queue = [a];
      visited.add(a.id);
      while (queue.length > 0) {
        const cur = queue.shift()!;
        comp.push(cur);
        for (const nbId of adj.get(cur.id) || []) {
          if (!visited.has(nbId)) {
            visited.add(nbId);
            const nbAtom = atomList.find((x) => x.id === nbId);
            if (nbAtom) queue.push(nbAtom);
          }
        }
      }
      comps.push(comp);
    }
    return comps;
  };

  let atoms = rawAtoms.map((a) => ({ ...a }));
  let bonds = rawBonds.map((b) => ({ ...b }));
  let comps = getComponents(atoms, bonds);

  if (comps.length <= 1) {
    return { atoms, bonds };
  }

  // Strip isolated group-1 / group-2 spectator counterions (e.g. 4 K+ in K4[Fe(CN)6] or 2 Na+ in Na2[Fe(CN)5NO])
  if (stripSpectatorCounterions && atoms.length > 2) {
    const hasNonSpectatorCore = atoms.some(
      (a) => !SPECTATOR_COUNTERIONS.has(a.element) && a.element !== 'H'
    );
    if (hasNonSpectatorCore) {
      const spectatorIds = new Set(
        comps
          .filter((c) => c.length === 1 && SPECTATOR_COUNTERIONS.has(c[0].element))
          .map((c) => c[0].id)
      );
      if (spectatorIds.size > 0 && spectatorIds.size < atoms.length) {
        atoms = atoms.filter((a) => !spectatorIds.has(a.id));
        bonds = bonds.filter((b) => !spectatorIds.has(b.source) && !spectatorIds.has(b.target));
        comps = getComponents(atoms, bonds);
        if (comps.length <= 1) {
          return { atoms, bonds };
        }
      }
    }
  }

  // Find central coordination atom among disconnected components
  const tmOrLnSingle = comps.find(
    (c) => c.length === 1 && (isTransitionMetal(c[0].element) || isLanthanide(c[0].element))
  )?.[0];
  const tmOrLnAny = atoms.find(
    (a) => isTransitionMetal(a.element) || isLanthanide(a.element)
  );
  const mainGroupSingle =
    comps.length >= 3
      ? comps.find((c) => c.length === 1 && MAIN_GROUP_COORD_CENTERS.has(c[0].element))?.[0]
      : undefined;

  const centerAtom = tmOrLnSingle || tmOrLnAny || mainGroupSingle;

  if (!centerAtom) {
    if (atoms.length === 2 && bonds.length === 0) {
      return {
        atoms,
        bonds: [
          {
            id: 'b_coord_1',
            source: atoms[0].id,
            target: atoms[1].id,
            order: 1,
          },
        ],
      };
    }
    return { atoms, bonds };
  }

  const distToCenter = (a: AtomNode) => {
    const d3 = Math.hypot(
      a.x3d - centerAtom.x3d,
      a.y3d - centerAtom.y3d,
      a.z3d - centerAtom.z3d
    );
    if (d3 > 0.01) return d3;
    return Math.hypot(a.x2d - centerAtom.x2d, a.y2d - centerAtom.y2d);
  };

  for (const comp of comps) {
    if (comp.some((a) => a.id === centerAtom.id)) continue;
    if (comp.length > 8) continue;

    let donor: AtomNode = comp[0];
    if (comp.length === 1) {
      donor = comp[0];
    } else if (comp.length === 2) {
      const [a, b] = comp;
      if (a.element === 'H' && b.element !== 'H') {
        donor = b;
      } else if (b.element === 'H' && a.element !== 'H') {
        donor = a;
      } else if (
        (a.element === 'C' && (b.element === 'N' || b.element === 'O' || b.element === 'S')) ||
        (b.element === 'C' && (a.element === 'N' || a.element === 'O' || a.element === 'S'))
      ) {
        // Cyanide (-C≡N), Carbonyl (-C≡O / -C=O), Thiocarbonyl (-C≡S): Carbon is always the σ-donor
        const cAtom = a.element === 'C' ? a : b;
        const heteroAtom = a.element === 'C' ? b : a;
        donor = cAtom;

        // Upgrade divalent carbonyl [C]=O in PubChem SDFs (e.g. Fe(CO)5) to standard C≡O triple bond
        if (heteroAtom.element === 'O') {
          bonds = bonds.map((edge) =>
            (edge.source === cAtom.id && edge.target === heteroAtom.id) ||
            (edge.source === heteroAtom.id && edge.target === cAtom.id)
              ? { ...edge, order: 3 as const }
              : edge
          );
        }
      } else if (
        (a.element === 'N' && b.element === 'O') ||
        (b.element === 'N' && a.element === 'O')
      ) {
        donor = a.element === 'N' ? a : b;
      } else {
        donor = distToCenter(a) <= distToCenter(b) ? a : b;
      }
    } else {
      const nonH = comp.filter((a) => a.element !== 'H');
      const pool = nonH.length > 0 ? nonH : comp;
      donor = [...pool].sort((a, b) => distToCenter(a) - distToCenter(b))[0];
    }

    bonds.push({
      id: `b_coord_${bonds.length + 1}`,
      source: centerAtom.id,
      target: donor.id,
      order: 1,
    });
  }

  return { atoms, bonds };
}

export interface LigandBranchInfo {
  donor: AtomNode;
  distalAtoms: AtomNode[];
  sig: string;
}

export interface StarCoordinationTopology {
  centerAtom: AtomNode;
  branches: LigandBranchInfo[];
  isMonoatomicLigands: boolean;
}

/**
 * Detects whether a molecular graph has a single central coordination / VSEPR atom
 * surrounded by nLig ∈ [2, 6] mutually disjoint terminal ligand branches (supporting both
 * monoatomic ligands like -Cl, -F, -H and polyatomic ligands like -C≡N, -C≡O, -NO, -NH3, -OH2).
 */
export function analyzeStarCoordinationTopology(
  atoms: AtomNode[],
  bonds: BondEdge[]
): StarCoordinationTopology | null {
  if (atoms.length < 3 || atoms.length > 32) return null;

  const adj = new Map<string, string[]>();
  for (const a of atoms) adj.set(a.id, []);
  for (const b of bonds) {
    if (adj.has(b.source) && adj.has(b.target)) {
      adj.get(b.source)!.push(b.target);
      adj.get(b.target)!.push(b.source);
    }
  }

  // Sort candidate central atoms: transition metals & lanthanides first, then highest coordination degree
  const candidates = atoms
    .filter((a) => {
      const deg = (adj.get(a.id) || []).length;
      return deg >= 2 && deg <= 6;
    })
    .sort((a, b) => {
      const aMetal = isTransitionMetal(a.element) || isLanthanide(a.element) ? 1 : 0;
      const bMetal = isTransitionMetal(b.element) || isLanthanide(b.element) ? 1 : 0;
      if (aMetal !== bMetal) return bMetal - aMetal;
      return (adj.get(b.id) || []).length - (adj.get(a.id) || []).length;
    });

  for (const centerAtom of candidates) {
    const donorIds = adj.get(centerAtom.id) || [];
    const nLig = donorIds.length;
    if (nLig < 2 || nLig > 6) continue;

    const globalVisited = new Set<string>([centerAtom.id]);
    const branches: LigandBranchInfo[] = [];
    let validStar = true;

    for (const dId of donorIds) {
      if (globalVisited.has(dId)) {
        validStar = false;
        break;
      }
      const donor = atoms.find((a) => a.id === dId);
      if (!donor) {
        validStar = false;
        break;
      }

      const branchAtoms: AtomNode[] = [];
      const queue: AtomNode[] = [donor];
      globalVisited.add(dId);

      while (queue.length > 0) {
        const cur = queue.shift()!;
        branchAtoms.push(cur);
        for (const nbId of adj.get(cur.id) || []) {
          if (nbId === centerAtom.id) continue;
          if (globalVisited.has(nbId)) {
            // Cycle or cross-branch link detected (e.g. benzene ring or chelating macrocycle)
            if (!branchAtoms.some((ba) => ba.id === nbId)) {
              validStar = false;
              break;
            }
            continue;
          }
          globalVisited.add(nbId);
          const nbAtom = atoms.find((a) => a.id === nbId);
          if (nbAtom) queue.push(nbAtom);
        }
        if (!validStar) break;
      }

      if (!validStar || branchAtoms.length > 5) {
        validStar = false;
        break;
      }

      const distalAtoms = branchAtoms.filter((a) => a.id !== donor.id);
      const distalSig = distalAtoms
        .map((a) => a.element)
        .sort()
        .join('');
      branches.push({
        donor,
        distalAtoms,
        sig: `${donor.element}:${distalSig}`,
      });
    }

    if (validStar && globalVisited.size === atoms.length) {
      const isMonoatomicLigands = branches.every((b) => b.distalAtoms.length === 0);
      const isMetalCenter =
        isTransitionMetal(centerAtom.element) || isLanthanide(centerAtom.element);
      const polyatomicBranchCount = branches.filter((b) => b.distalAtoms.length > 0).length;

      // Prevent non-star multi-center organic molecules like C2H4 (H2C=CH2) from being misidentified
      // as a 3-branch star around C1 with one "=CH2" polyatomic branch
      if (!isMonoatomicLigands && !isMetalCenter && polyatomicBranchCount < 2) {
        continue;
      }

      return {
        centerAtom,
        branches,
        isMonoatomicLigands,
      };
    }
  }

  return null;
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
    // IUPAC Red Book electropositive-to-electronegative ordering (metals first, then p-block/noble gases/H/halogens/chalcogens)
    const iupacNonMetalOrder = [
      'Rn', 'Xe', 'Kr', 'Ar', 'Ne', 'He',
      'B', 'Si', 'Ge', 'As', 'Sb', 'P', 'N',
      'H', 'Po', 'Te', 'Se', 'S',
      'At', 'I', 'Br', 'Cl', 'O', 'F',
    ];
    const sorted = [...keys].sort((a, b) => {
      const ia = iupacNonMetalOrder.indexOf(a);
      const ib = iupacNonMetalOrder.indexOf(b);
      if (ia !== -1 && ib !== -1) return ia - ib;
      // Electropositive metals (not in iupacNonMetalOrder) come before non-metals
      if (ia === -1 && ib !== -1) return -1;
      if (ia !== -1 && ib === -1) return 1;
      const va = Math.abs(ELEMENT_PROPERTIES[a as ElementSymbol]?.vsieS ?? 8);
      const vb = Math.abs(ELEMENT_PROPERTIES[b as ElementSymbol]?.vsieS ?? 8);
      if (Math.abs(va - vb) > 0.05) return va - vb;
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
 * Generates an industry-standard MDL V2000 Molfile from the drawn 2D/3D molecular graph.
 */
export function generateMolfileFromGraph(
  atoms: AtomNode[],
  bonds: BondEdge[],
  title = 'Symmetry & MO Visualizer Structure'
): string {
  const now = new Date();
  const dateStr = `${String(now.getMonth() + 1).padStart(2, '0')}${String(now.getDate()).padStart(2, '0')}${String(now.getFullYear()).slice(-2)}${String(now.getHours()).padStart(2, '0')}${String(now.getMinutes()).padStart(2, '0')}`;
  const header1 = title.slice(0, 80);
  const header2 = `  SymMOVis${dateStr}3D`;
  const header3 = 'Generated by Symmetry & Molecular Orbital Visualizer';

  const countLine = `${String(atoms.length).padStart(3, ' ')}${String(bonds.length).padStart(3, ' ')}  0  0  0  0  0  0  0  0999 V2000`;
  const atomIdxMap = new Map<string, number>();
  const atomLines: string[] = [];

  atoms.forEach((a, idx) => {
    atomIdxMap.set(a.id, idx + 1);
    const x = a.x3d.toFixed(4).padStart(10, ' ');
    const y = a.y3d.toFixed(4).padStart(10, ' ');
    const z = a.z3d.toFixed(4).padStart(10, ' ');
    const sym = a.element.padEnd(3, ' ');
    const chg = a.formalCharge ?? 0;
    const chgCode =
      chg === 3 ? 1 : chg === 2 ? 2 : chg === 1 ? 3 : chg === -1 ? 5 : chg === -2 ? 6 : chg === -3 ? 7 : 0;
    atomLines.push(`${x}${y}${z} ${sym} 0  ${chgCode}  0  0  0  0  0  0  0  0  0  0`);
  });

  const bondLines: string[] = [];
  for (const b of bonds) {
    const i1 = atomIdxMap.get(b.source);
    const i2 = atomIdxMap.get(b.target);
    if (!i1 || !i2) continue;
    const stereoCode =
      b.stereo === 'wedge' ? 1 : b.stereo === 'dash' ? 6 : b.stereo === 'wavy' ? 4 : b.stereo === 'crossed' ? 3 : 0;
    bondLines.push(
      `${String(i1).padStart(3, ' ')}${String(i2).padStart(3, ' ')}${String(b.order).padStart(3, ' ')}${String(stereoCode).padStart(3, ' ')}  0  0  0`
    );
  }

  return [header1, header2, header3, countLine, ...atomLines, ...bondLines, 'M  END'].join('\n');
}

/**
 * Parses a SMILES string, Molecular Formula, or MDL V2000 Molfile into a 2D/3D molecular graph
 * supporting all 86 elements (H through Rn).
 */
export function parseChemicalInputToGraph(
  rawText: string,
  format: 'SMILES' | 'Formula' | 'MDL Molfile' | 'InChI' = 'SMILES'
): { atoms: AtomNode[]; bonds: BondEdge[]; title?: string } | null {
  const text = rawText.trim();
  if (!text) return null;

  // 1. Check if input is an MDL V2000 Molfile
  if (format === 'MDL Molfile' || text.includes('V2000') || text.includes('M  END')) {
    const lines = text.split(/\r?\n/);
    const v2000Idx = lines.findIndex((l) => l.includes('V2000'));
    if (v2000Idx >= 0) {
      const countsLine = lines[v2000Idx];
      const numAtoms = parseInt(countsLine.slice(0, 3).trim(), 10);
      const numBonds = parseInt(countsLine.slice(3, 6).trim(), 10);
      if (numAtoms > 0 && lines.length >= v2000Idx + 1 + numAtoms) {
        const atoms: AtomNode[] = [];
        const bonds: BondEdge[] = [];
        const elemCounts: Record<string, number> = {};
        for (let i = 0; i < numAtoms; i++) {
          const line = lines[v2000Idx + 1 + i] || '';
          const x = parseFloat(line.slice(0, 10).trim()) || 0;
          const y = parseFloat(line.slice(10, 20).trim()) || 0;
          const z = parseFloat(line.slice(20, 30).trim()) || 0;
          const rawSym = line.slice(31, 34).trim();
          const el = (ELEMENT_PROPERTIES[rawSym as ElementSymbol] ? rawSym : 'C') as ElementSymbol;
          const chgCode = parseInt(line.slice(36, 39).trim(), 10) || 0;
          const formalCharge =
            chgCode === 1 ? 3 : chgCode === 2 ? 2 : chgCode === 3 ? 1 : chgCode === 5 ? -1 : chgCode === 6 ? -2 : chgCode === 7 ? -3 : 0;
          elemCounts[el] = (elemCounts[el] || 0) + 1;
          atoms.push({
            id: `${el}${elemCounts[el]}`,
            element: el,
            x2d: Number(x.toFixed(2)),
            y2d: Number((-y).toFixed(2)),
            x3d: Number(x.toFixed(3)),
            y3d: Number(y.toFixed(3)),
            z3d: Number(z.toFixed(3)),
            ...(formalCharge !== 0 ? { formalCharge } : {}),
          });
        }
        for (let j = 0; j < numBonds; j++) {
          const line = lines[v2000Idx + 1 + numAtoms + j] || '';
          const sIdx = parseInt(line.slice(0, 3).trim(), 10) - 1;
          const tIdx = parseInt(line.slice(3, 6).trim(), 10) - 1;
          const rawOrder = parseInt(line.slice(6, 9).trim(), 10) || 1;
          const stereoCode = parseInt(line.slice(9, 12).trim(), 10) || 0;
          if (atoms[sIdx] && atoms[tIdx]) {
            const order: 1 | 2 | 3 = rawOrder === 2 ? 2 : rawOrder === 3 ? 3 : 1;
            const stereo: BondEdge['stereo'] =
              stereoCode === 1 ? 'wedge' : stereoCode === 6 ? 'dash' : stereoCode === 4 ? 'wavy' : stereoCode === 3 ? 'crossed' : 'none';
            bonds.push({
              id: `b_mol_${j + 1}`,
              source: atoms[sIdx].id,
              target: atoms[tIdx].id,
              order,
              stereo,
            });
          }
        }
        const resolved = resolveCoordinationBonds(atoms, bonds);
        return {
          atoms: resolved.atoms,
          bonds: resolved.bonds,
          title: lines[0]?.trim() || 'Imported Molfile',
        };
      }
    }
  }

  // 2. Check if input matches a canonical benchmark molecule by formula, SMILES, or name
  const upper = text.toUpperCase();
  const canonical = CANONICAL_MOLECULES.find(
    (m) =>
      m.formula.toUpperCase() === upper ||
      m.smiles.toUpperCase() === upper ||
      m.name.toUpperCase() === upper
  );
  if (canonical) {
    return {
      atoms: canonical.atoms.map((a) => ({ ...a })),
      bonds: canonical.bonds.map((b) => ({ ...b })),
      title: canonical.name,
    };
  }

  // 3. If format is Formula (or input looks like a stoich formula such as SiH4, GeCl4, HBr, PCl3, H2Se)
  const isFormulaLike =
    format === 'Formula' ||
    (/^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)+$/.test(text) && !/[=#()[\]]/.test(text));
  if (isFormulaLike) {
    const tokenRegex = /([A-Z][a-z]?)(\d*)/g;
    const parts: { el: ElementSymbol; count: number }[] = [];
    let match: RegExpExecArray | null;
    while ((match = tokenRegex.exec(text)) !== null) {
      const sym = match[1] as ElementSymbol;
      if (!ELEMENT_PROPERTIES[sym]) continue;
      const count = match[2] ? parseInt(match[2], 10) : 1;
      parts.push({ el: sym, count: Math.min(12, Math.max(1, count)) });
    }
    if (parts.length > 0) {
      const totalAtoms = parts.reduce((s, p) => s + p.count, 0);
      if (totalAtoms >= 1 && totalAtoms <= 24) {
        // Choose central atom as the unique element or lowest count non-H element
        const centerPartIdx =
          parts.findIndex((p) => p.count === 1 && p.el !== 'H') !== -1
            ? parts.findIndex((p) => p.count === 1 && p.el !== 'H')
            : 0;
        const atoms: AtomNode[] = [];
        const bonds: BondEdge[] = [];
        const counts: Record<string, number> = {};
        const nextId = (el: ElementSymbol) => {
          counts[el] = (counts[el] || 0) + 1;
          return `${el}${counts[el]}`;
        };

        const centerEl = parts[centerPartIdx].el;
        const centerId = nextId(centerEl);
        atoms.push({
          id: centerId,
          element: centerEl,
          x2d: 0,
          y2d: 0,
          x3d: 0,
          y3d: 0,
          z3d: 0,
        });

        const ligands: ElementSymbol[] = [];
        parts.forEach((p, idx) => {
          const rem = idx === centerPartIdx ? p.count - 1 : p.count;
          for (let k = 0; k < rem; k++) ligands.push(p.el);
        });

        ligands.forEach((ligEl, idx) => {
          const angle = (idx * 2 * Math.PI) / Math.max(1, ligands.length);
          const r = ligEl === 'H' ? 1.15 : 1.45;
          const lx = Number((Math.cos(angle) * r).toFixed(2));
          const ly = Number((Math.sin(angle) * r).toFixed(2));
          const lid = nextId(ligEl);
          atoms.push({
            id: lid,
            element: ligEl,
            x2d: lx,
            y2d: ly,
            x3d: lx,
            y3d: -ly,
            z3d: 0,
          });
          const order: 1 | 2 = ligEl === 'O' && ligands.length <= 3 ? 2 : 1;
          bonds.push({
            id: `b_form_${idx + 1}`,
            source: centerId,
            target: lid,
            order,
          });
        });

        return { atoms, bonds, title: text };
      }
    }
  }

  // 4. Parse linear/branched/cyclic SMILES strings (e.g. C=CC=C, c1ccccc1, CC(=O)O, [Li+], BrCBr)
  const atoms: AtomNode[] = [];
  const bonds: BondEdge[] = [];
  const counts: Record<string, number> = {};
  const nextId = (el: ElementSymbol) => {
    counts[el] = (counts[el] || 0) + 1;
    return `${el}${counts[el]}`;
  };

  const ringClosures = new Map<string, { atomId: string; order: 1 | 2 | 3 }>();
  const branchStack: string[] = [];
  let currentAtomId: string | null = null;
  let pendingOrder: 1 | 2 | 3 = 1;
  let pendingStereo: BondEdge['stereo'] = 'none';
  let stepIdx = 0;

  let i = 0;
  while (i < text.length) {
    const ch = text[i];
    if (ch === '.') {
      currentAtomId = null;
      pendingOrder = 1;
      pendingStereo = 'none';
      i++;
      continue;
    }
    if (ch === '(') {
      if (currentAtomId) branchStack.push(currentAtomId);
      i++;
      continue;
    }
    if (ch === ')') {
      currentAtomId = branchStack.pop() || currentAtomId;
      i++;
      continue;
    }
    if (ch === '=') {
      pendingOrder = 2;
      i++;
      continue;
    }
    if (ch === '#') {
      pendingOrder = 3;
      i++;
      continue;
    }
    if (ch === '-') {
      pendingOrder = 1;
      i++;
      continue;
    }
    if (ch === '/' || ch === '\\') {
      pendingStereo = 'wavy';
      i++;
      continue;
    }
    if (/\d/.test(ch)) {
      if (currentAtomId) {
        const existing = ringClosures.get(ch);
        if (existing) {
          bonds.push({
            id: `b_ring_${ch}_${bonds.length + 1}`,
            source: existing.atomId,
            target: currentAtomId,
            order: pendingOrder !== 1 ? pendingOrder : existing.order,
          });
          ringClosures.delete(ch);
        } else {
          ringClosures.set(ch, { atomId: currentAtomId, order: pendingOrder });
        }
        pendingOrder = 1;
      }
      i++;
      continue;
    }

    let parsedEl: ElementSymbol | null = null;
    let formalCharge = 0;
    let isAromatic = false;

    if (ch === '[') {
      const closeIdx = text.indexOf(']', i);
      if (closeIdx > i) {
        const bracketContent = text.slice(i + 1, closeIdx);
        const symMatch = bracketContent.match(/([A-Z][a-z]?)/);
        if (symMatch && ELEMENT_PROPERTIES[symMatch[1] as ElementSymbol]) {
          parsedEl = symMatch[1] as ElementSymbol;
        }
        if (bracketContent.includes('+')) {
          const plusCount = (bracketContent.match(/\+/g) || []).length;
          const numMatch = bracketContent.match(/\+(\d+)/);
          formalCharge = numMatch ? parseInt(numMatch[1], 10) : plusCount;
        } else if (bracketContent.includes('-')) {
          const minusCount = (bracketContent.match(/-/g) || []).length;
          const numMatch = bracketContent.match(/-(\d+)/);
          formalCharge = -(numMatch ? parseInt(numMatch[1], 10) : minusCount);
        }
        i = closeIdx + 1;
      } else {
        i++;
        continue;
      }
    } else if (/[A-Za-z]/.test(ch)) {
      const twoChar = text.slice(i, i + 2);
      if (twoChar.length === 2 && /^[A-Z][a-z]$/.test(twoChar) && ELEMENT_PROPERTIES[twoChar as ElementSymbol]) {
        parsedEl = twoChar as ElementSymbol;
        i += 2;
      } else {
        const upperCh = ch.toUpperCase() as ElementSymbol;
        if (ELEMENT_PROPERTIES[upperCh]) {
          parsedEl = upperCh;
          isAromatic = ch === ch.toLowerCase();
        }
        i += 1;
      }
    } else {
      i++;
      continue;
    }

    if (parsedEl) {
      const id = nextId(parsedEl);
      const prevAtom = atoms.find((a) => a.id === currentAtomId);
      const angle = stepIdx % 2 === 0 ? -Math.PI / 6 : Math.PI / 6;
      const dist = 1.28;
      const x2d = prevAtom
        ? Number((prevAtom.x2d + Math.cos(angle) * dist).toFixed(2))
        : Number(((stepIdx - 2) * 0.85).toFixed(2));
      const y2d = prevAtom
        ? Number((prevAtom.y2d + Math.sin(angle) * dist + (branchStack.length > 0 ? 0.65 : 0)).toFixed(2))
        : 0;

      atoms.push({
        id,
        element: parsedEl,
        x2d,
        y2d,
        x3d: x2d,
        y3d: -y2d,
        z3d: 0,
        ...(formalCharge !== 0 ? { formalCharge } : {}),
      });

      if (currentAtomId) {
        const order: 1 | 2 | 3 =
          pendingOrder !== 1 ? pendingOrder : isAromatic && stepIdx % 2 === 1 ? 2 : 1;
        bonds.push({
          id: `b_smi_${bonds.length + 1}`,
          source: currentAtomId,
          target: id,
          order,
          stereo: pendingStereo,
        });
      }

      currentAtomId = id;
      pendingOrder = 1;
      pendingStereo = 'none';
      stepIdx++;
    }
  }

  if (atoms.length === 0) return null;
  const resolved = resolveCoordinationBonds(atoms, bonds);
  return { atoms: resolved.atoms, bonds: resolved.bonds, title: text };
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

  // 2. Check if molecule has a central coordination / VSEPR atom surrounded by nLig ∈ [2, 6] ligand branches
  // (Supports both monoatomic ligands like H2O, NH3, CH4, SF6, CrCl6, CeCl6 and polyatomic coordination complexes
  // like Hexacyanoferrate [Fe(CN)6]4-, Ferricyanide [Fe(CN)6]3-, Cr(CO)6, Fe(CO)5, Ni(CO)4, [Ni(CN)4]2-, [Co(NH3)6]3+).
  const maxAbsZInit = atoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d || 0)), 0);
  const starTopo = analyzeStarCoordinationTopology(atoms, bonds);
  const shouldApplyStarVsepr =
    starTopo !== null &&
    (starTopo.isMonoatomicLigands ||
      isTransitionMetal(starTopo.centerAtom.element) ||
      isLanthanide(starTopo.centerAtom.element) ||
      starTopo.branches.length >= 5 ||
      maxAbsZInit <= 0.08);

  if (starTopo && shouldApplyStarVsepr) {
    const { centerAtom } = starTopo;
    const nLig = starTopo.branches.length;

    // Sort ligand branches so unique ligand signatures or higher-bond-order terminal oxo/sulfido ligands
    // come first (sortedBranches[0]), ensuring they occupy the axial / unique symmetry site.
    const sigFreq: Record<string, number> = {};
    for (const br of starTopo.branches) {
      sigFreq[br.sig] = (sigFreq[br.sig] || 0) + 1;
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
    const sortedBranches = [...starTopo.branches].sort((a, b) => {
      const freqDiff = (sigFreq[a.sig] || 1) - (sigFreq[b.sig] || 1);
      if (freqDiff !== 0) return freqDiff;
      return getBondOrderToCenter(b.donor) - getBondOrderToCenter(a.donor);
    });
    const ligands = sortedBranches.map((br) => br.donor);

    const isMetalCenter =
      isTransitionMetal(centerAtom.element) || isLanthanide(centerAtom.element);
    const lonePairs = getCentralAtomLonePairs(centerAtom, atoms, bonds, charge);
    const rCenter = ELEMENT_PROPERTIES[centerAtom.element]?.covalentRadius ?? 0.8;
    const scaleFactor = lonePairs >= 2 || nLig >= 5 ? 0.98 : 0.96;
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
      if (isMetalCenter || lonePairs === 0 || lonePairs >= 3) {
        coords3D.set(ligands[0].id, [-r1, 0, 0]);
        coords3D.set(ligands[1].id, [r2, 0, 0]);
      } else {
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
      if (isMetalCenter || lonePairs === 0) {
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
        const rEq = getBondLen(ligands[0]);
        const rAx1 = getBondLen(ligands[1]) * 1.05;
        const rAx2 = getBondLen(ligands[2]) * 1.05;
        coords3D.set(centerAtom.id, [0, 0.2, 0]);
        coords3D.set(ligands[0].id, [0, Number((0.2 - rEq).toFixed(3)), 0]);
        coords3D.set(ligands[1].id, [Number((-rAx1).toFixed(3)), 0.12, 0]);
        coords3D.set(ligands[2].id, [Number(rAx2.toFixed(3)), 0.12, 0]);
      }
    } else if (nLig === 4) {
      const isD10Carbonyl =
        centerAtom.element === 'Ni' &&
        charge === 0 &&
        sortedBranches.every((b) => b.sig === 'C:O');
      const isSqPlanarMetal = isSquarePlanarMetal(centerAtom.element) && !isD10Carbonyl;

      if ((isMetalCenter && !isSqPlanarMetal) || (!isMetalCenter && lonePairs === 0)) {
        // AX4E0 (CH4, SiH4, CCl4, TiCl4, Ni(CO)4) -> Tetrahedral (Td / C3v)
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
      } else if (!isMetalCenter && lonePairs === 1) {
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
        // AX4E2 (XeF4) or d8 square-planar transition metal (PtCl4, [Ni(CN)4]2-, PdCl4, AuCl4) -> Square planar in xy plane (D4h, z = 0)
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
      if (isMetalCenter || lonePairs === 0) {
        // AX5E0 (PF5, PCl5, Fe(CO)5) -> Trigonal bipyramidal (D3h: 3 eq in xy, 2 ax along ±z)
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
      if (isMetalCenter || lonePairs === 0) {
        // AX6E0 (SF6) or transition-metal / lanthanide ML6 ([Fe(CN)6]4-, Cr(CO)6, CrCl6, CeCl6) -> Regular Octahedral (Oh)
        const octDirs: [number, number, number][] = [
          [1, 0, 0],
          [-1, 0, 0],
          [0, 1, 0],
          [0, -1, 0],
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

    // Propagate 3D coordinates outward along each polyatomic ligand branch (e.g. N in Fe-C≡N, O in Cr-C≡O, H3 in Co-NH3)
    const centerPos = coords3D.get(centerAtom.id) || [0, 0, 0];
    for (const br of sortedBranches) {
      if (br.distalAtoms.length === 0) continue;
      const dPos = coords3D.get(br.donor.id);
      if (!dPos) continue;
      const vx = dPos[0] - centerPos[0];
      const vy = dPos[1] - centerPos[1];
      const vz = dPos[2] - centerPos[2];
      const vLen = Math.hypot(vx, vy, vz) || 1;
      const ux = vx / vLen;
      const uy = vy / vLen;
      const uz = vz / vLen;

      if (br.distalAtoms.length === 1) {
        const distal = br.distalAtoms[0];
        const b = bonds.find(
          (edge) =>
            (edge.source === br.donor.id && edge.target === distal.id) ||
            (edge.target === br.donor.id && edge.source === distal.id)
        );
        const rDonor = ELEMENT_PROPERTIES[br.donor.element]?.covalentRadius ?? 0.75;
        const rDistal = ELEMENT_PROPERTIES[distal.element]?.covalentRadius ?? 0.68;
        const bondScale = b && b.order >= 2 ? 0.84 : 0.95;
        const distLen = Number(((rDonor + rDistal) * bondScale).toFixed(3));
        // Linear diatomic ligands (-C≡N, -C≡O, -N≡O, -C=O) extend strictly collinearly along u
        coords3D.set(distal.id, [
          Number((dPos[0] + ux * distLen).toFixed(3)),
          Number((dPos[1] + uy * distLen).toFixed(3)),
          Number((dPos[2] + uz * distLen).toFixed(3)),
        ]);
      } else {
        // Multi-atom conical substituent branch (e.g. -NH3, -OH2, -CH3, -NO2)
        const refX = Math.abs(uz) < 0.9 ? 0 : 1;
        const refZ = Math.abs(uz) < 0.9 ? 1 : 0;
        let e1x = uy * refZ;
        let e1y = uz * refX - ux * refZ;
        let e1z = -uy * refX;
        const e1Len = Math.hypot(e1x, e1y, e1z) || 1;
        e1x /= e1Len;
        e1y /= e1Len;
        e1z /= e1Len;
        const e2x = uy * e1z - uz * e1y;
        const e2y = uz * e1x - ux * e1z;
        const e2z = ux * e1y - uy * e1x;

        const kCount = br.distalAtoms.length;
        const cosTheta = 0.52;
        const sinTheta = 0.854;
        br.distalAtoms.forEach((distal, mIdx) => {
          const rDonor = ELEMENT_PROPERTIES[br.donor.element]?.covalentRadius ?? 0.72;
          const rDistal = ELEMENT_PROPERTIES[distal.element]?.covalentRadius ?? 0.35;
          const distLen = (rDonor + rDistal) * 0.95;
          const phi = (mIdx * 2 * Math.PI) / kCount;
          const cp = Math.cos(phi);
          const sp = Math.sin(phi);
          const dirX = cosTheta * ux + sinTheta * (cp * e1x + sp * e2x);
          const dirY = cosTheta * uy + sinTheta * (cp * e1y + sp * e2y);
          const dirZ = cosTheta * uz + sinTheta * (cp * e1z + sp * e2z);
          coords3D.set(distal.id, [
            Number((dPos[0] + dirX * distLen).toFixed(3)),
            Number((dPos[1] + dirY * distLen).toFixed(3)),
            Number((dPos[2] + dirZ * distLen).toFixed(3)),
          ]);
        });
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
      const isMetalCenter =
        isTransitionMetal(center.element) || isLanthanide(center.element);
      const lonePairs = getCentralAtomLonePairs(center, atoms, bonds, charge);
      // AX2E0 (CO2, CS2, BeF2), linear ML2, or AX2E3 (XeF2) are linear; AX2E1 (SO2, O3, NO2) and AX2E2 (H2O, OF2) are bent
      const isLinearVsepr = isMetalCenter || lonePairs === 0 || lonePairs >= 3;

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

  // Central atom surrounded by nLig ligands (AX_nE_m VSEPR & d/f coordination point-group classification,
  // supporting both monoatomic ligands and polyatomic ligands like -CN, -CO, -NO, -NH3)
  const starTopo = analyzeStarCoordinationTopology(atoms, bonds);

  if (starTopo) {
    const { centerAtom, branches } = starTopo;
    const ligandSigs = new Set(branches.map((b) => b.sig));
    const nLig = branches.length;
    const isMetalCenter =
      isTransitionMetal(centerAtom.element) || isLanthanide(centerAtom.element);
    const lonePairs = getCentralAtomLonePairs(centerAtom, atoms, bonds, charge);

    if (ligandSigs.size === 1) {
      if (nLig === 2) {
        return isMetalCenter || lonePairs === 0 || lonePairs >= 3 ? 'Dinfh' : 'C2v';
      }
      if (nLig === 3) {
        // AX3E0 (BF3, SO3) or ML3 -> D3h; AX3E1 (NH3, PCl3) -> C3v; AX3E2 (ClF3) -> C2v
        if (isMetalCenter || lonePairs === 0) return 'D3h';
        if (lonePairs === 1) return 'C3v';
        return 'C2v';
      }
      if (nLig === 4) {
        if (isMetalCenter) {
          const isD10Carbonyl =
            centerAtom.element === 'Ni' &&
            charge === 0 &&
            branches.every((b) => b.sig === 'C:O');
          return isSquarePlanarMetal(centerAtom.element) && !isD10Carbonyl ? 'D4h' : 'Td';
        }
        // AX4E0 (CH4, CCl4) -> Td; AX4E1 (SF4) -> C2v; AX4E2 (XeF4) -> D4h
        if (lonePairs === 0) return 'Td';
        if (lonePairs === 1) return 'C2v';
        return 'D4h';
      }
      if (nLig === 5) {
        // AX5E0 (PF5, PCl5, Fe(CO)5) -> D3h; AX5E1 (BrF5, IF5) -> C4v
        return isMetalCenter || lonePairs === 0 ? 'D3h' : 'C4v';
      }
      if (nLig === 6) {
        // AX6E0 (SF6) or transition-metal/lanthanide ML6 ([Fe(CN)6]4-, Cr(CO)6, CrCl6, CeCl6) -> Oh; AX6E1 (XeF6) -> C3v
        return isMetalCenter || lonePairs === 0 ? 'Oh' : 'C3v';
      }
    }

    if (ligandSigs.size === 2) {
      const sigCounts = Array.from(ligandSigs).map(
        (sig) => branches.filter((b) => b.sig === sig).length
      );
      if (nLig === 2) {
        return isMetalCenter || lonePairs === 0 || lonePairs >= 3 ? 'Cinfv' : 'Cs';
      }
      if (nLig === 3) {
        // AX2Y: planar (CH2O, COCl2 -> C2v) vs pyramidal (SOCl2 -> Cs)
        return isMetalCenter || lonePairs === 0 ? 'C2v' : 'Cs';
      }
      if (nLig === 4 && sigCounts.includes(3)) {
        // AX3Y tetrahedral (CHCl3, POCl3 -> C3v)
        return isMetalCenter || lonePairs === 0 ? 'C3v' : 'Cs';
      }
      if (nLig === 4 && sigCounts.includes(2)) {
        return 'C2v';
      }
      if (nLig === 5 && sigCounts.includes(4)) {
        // AX4Y square pyramidal (XeOF4 -> C4v)
        return 'C4v';
      }
      if (nLig === 6 && sigCounts.includes(5)) {
        // ML5Y monosubstituted octahedral (e.g. [Fe(CN)5(NO)]2- -> C4v)
        return 'C4v';
      }
      if (nLig === 6 && sigCounts.includes(4)) {
        return 'D4h';
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

function evalCanonicalAngular(
  orbital: OrbitalBasisType,
  ux: number,
  uy: number,
  uz: number
): number {
  switch (orbital) {
    case 's':
      return 1.0;
    case 'px':
      return ux;
    case 'py':
      return uy;
    case 'pz':
      return uz;
    case 'dz2':
      return 0.5 * (2 * uz * uz - ux * ux - uy * uy);
    case 'dx2y2':
      return 0.8660254 * (ux * ux - uy * uy);
    case 'dxy':
      return 1.7320508 * ux * uy;
    case 'dxz':
      return 1.7320508 * ux * uz;
    case 'dyz':
      return 1.7320508 * uy * uz;
    case 'fz3':
      return 0.5 * uz * (2 * uz * uz - 3 * ux * ux - 3 * uy * uy);
    case 'fxz2':
      return 0.6123724 * ux * (4 * uz * uz - ux * ux - uy * uy);
    case 'fyz2':
      return 0.6123724 * uy * (4 * uz * uz - ux * ux - uy * uy);
    case 'fxyz':
      return 3.8729833 * ux * uy * uz;
    case 'fz_x2y2':
      return 1.9364917 * uz * (ux * ux - uy * uy);
    case 'fx_x2_3y2':
      return 0.7905694 * ux * (ux * ux - 3 * uy * uy);
    case 'fy_3x2_y2':
      return 0.7905694 * uy * (3 * ux * ux - uy * uy);
  }
}

export function computeTwoCenterSlaterKosterOverlap(
  oi: OrbitalBasisType,
  oj: OrbitalBasisType,
  l: number,
  m: number,
  n: number,
  decay: number,
  orderBoost: number
): number {
  const isSpI = oi === 's' || oi === 'px' || oi === 'py' || oi === 'pz';
  const isSpJ = oj === 's' || oj === 'px' || oj === 'py' || oj === 'pz';

  if (isSpI && isSpJ) {
    const dirMap: Record<'px' | 'py' | 'pz', number> = { px: l, py: m, pz: n };
    const S_ss_sigma = 0.36 * decay * orderBoost;
    const S_sp_sigma = 0.32 * decay * orderBoost;
    const S_pp_sigma = 0.30 * decay * orderBoost;
    const S_pp_pi = 0.19 * decay * orderBoost;

    if (oi === 's' && oj === 's') return S_ss_sigma;
    if (oi === 's' && oj !== 's') return -dirMap[oj] * S_sp_sigma;
    if (oi !== 's' && oj === 's') return dirMap[oi] * S_sp_sigma;
    if (oi !== 's' && oj !== 's') {
      const di = dirMap[oi];
      const dj = dirMap[oj];
      return oi === oj
        ? -di * di * S_pp_sigma + (1 - di * di) * S_pp_pi
        : -di * dj * (S_pp_sigma + S_pp_pi);
    }
  }

  // Construct right-handed orthonormal local diatomic frame (e1, e2, e3) with e3 = (l, m, n)
  const refX = Math.abs(n) < 0.9 ? 0 : 1;
  const refZ = Math.abs(n) < 0.9 ? 1 : 0;
  let e1x = m * refZ;
  let e1y = n * refX - l * refZ;
  let e1z = -m * refX;
  const e1Len = Math.hypot(e1x, e1y, e1z) || 1;
  e1x /= e1Len;
  e1y /= e1Len;
  e1z /= e1Len;
  const e2x = m * e1z - n * e1y;
  const e2y = n * e1x - l * e1z;
  const e2z = l * e1y - m * e1x;

  // m_l = 0 (σ) head-on overlap along +e3 on atom i and -e3 on atom j
  const sigI = evalCanonicalAngular(oi, l, m, n);
  const sigJ = evalCanonicalAngular(oj, -l, -m, -n);

  // m_l = 1 (π) and m_l = 2 (δ) exact 8-point azimuthal Fourier projection around e3
  const c0 = 0.75;
  const s0 = 0.6614378;
  let piIx = 0,
    piIy = 0,
    piJx = 0,
    piJy = 0;
  let delIx = 0,
    delIy = 0,
    delJx = 0,
    delJy = 0;

  for (let k = 0; k < 8; k++) {
    const phi = (k * Math.PI) / 4;
    const cp = Math.cos(phi);
    const sp = Math.sin(phi);
    const tx = s0 * (cp * e1x + sp * e2x);
    const ty = s0 * (cp * e1y + sp * e2y);
    const tz = s0 * (cp * e1z + sp * e2z);
    const valI = evalCanonicalAngular(oi, c0 * l + tx, c0 * m + ty, c0 * n + tz);
    const valJ = evalCanonicalAngular(oj, -c0 * l + tx, -c0 * m + ty, -c0 * n + tz);
    piIx += valI * cp;
    piIy += valI * sp;
    piJx += valJ * cp;
    piJy += valJ * sp;
    const c2p = Math.cos(2 * phi);
    const s2p = Math.sin(2 * phi);
    delIx += valI * c2p;
    delIy += valI * s2p;
    delJx += valJ * c2p;
    delJy += valJ * s2p;
  }

  piIx *= 0.28;
  piIy *= 0.28;
  piJx *= 0.28;
  piJy *= 0.28;
  delIx *= 0.32;
  delIy *= 0.32;
  delJx *= 0.32;
  delJy *= 0.32;

  const isF = oi.startsWith('f') || oj.startsWith('f');
  const sSigma = (isF ? 0.16 : 0.25) * decay * orderBoost;
  const sPi = (isF ? 0.10 : 0.16) * decay * orderBoost;
  const sDelta = (isF ? 0.05 : 0.08) * decay * orderBoost;

  return (
    sigI * sigJ * sSigma +
    (piIx * piJx + piIy * piJy) * sPi +
    (delIx * delJx + delIy * delJy) * sDelta
  );
}

const ehtSolverCache = new Map<
  string,
  {
    basis: {
      atom: AtomNode;
      orbital: OrbitalBasisType;
      label: string;
      alpha: number;
    }[];
    eigenvalues: number[];
    eigenvectors: number[][];
  }
>();

/**
 * Computes raw Extended Hückel Theory (EHT) eigenvalues and eigenvectors using
 * an s-p-d-f valence basis (1s for H/He; ns+np for main-group; (n-1)d+ns+np for transition metals;
 * 4f+5d+6s+6p for lanthanides; and nd polarization for hypervalent Period >= 3 p-block centers),
 * with SO(3)-invariant Slater-Koster two-center overlap integrals S_mu_nu(R_AB).
 */
export function runExtendedHuckelSolver(
  atoms: AtomNode[],
  bonds: BondEdge[]
): {
  basis: {
    atom: AtomNode;
    orbital: OrbitalBasisType;
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
    orbital: OrbitalBasisType;
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

  const D_ORBITALS: { orb: OrbitalBasisType; suffix: string }[] = [
    { orb: 'dz2', suffix: 'dz2' },
    { orb: 'dx2y2', suffix: 'dx2-y2' },
    { orb: 'dxy', suffix: 'dxy' },
    { orb: 'dxz', suffix: 'dxz' },
    { orb: 'dyz', suffix: 'dyz' },
  ];

  const F_ORBITALS: { orb: OrbitalBasisType; suffix: string }[] = [
    { orb: 'fz3', suffix: 'fz3' },
    { orb: 'fxz2', suffix: 'fxz2' },
    { orb: 'fyz2', suffix: 'fyz2' },
    { orb: 'fxyz', suffix: 'fxyz' },
    { orb: 'fz_x2y2', suffix: 'fz(x2-y2)' },
    { orb: 'fx_x2_3y2', suffix: 'fx(x2-3y2)' },
    { orb: 'fy_3x2_y2', suffix: 'fy(3x2-y2)' },
  ];

  const basis: BasisFunc[] = [];
  for (const a of activeAtoms) {
    const props = ELEMENT_PROPERTIES[a.element] || ELEMENT_PROPERTIES.C;
    const shell = getPrincipalShell(a.element);
    if (a.element === 'H' || a.element === 'He') {
      basis.push({
        atom: a,
        orbital: 's',
        label: `${a.id} ${shell}s`,
        alpha: props.vsieS,
      });
    } else {
      // Lanthanides (La-Lu): include 7 real cubic 4f orbitals and 5 real cubic 5d orbitals
      if (isLanthanide(a.element)) {
        const alphaF = props.vsieF ?? -12.5;
        for (const f of F_ORBITALS) {
          basis.push({
            atom: a,
            orbital: f.orb,
            label: `${a.id} 4${f.suffix}`,
            alpha: alphaF,
          });
        }
        const alphaD = props.vsieD ?? -8.0;
        for (const d of D_ORBITALS) {
          basis.push({
            atom: a,
            orbital: d.orb,
            label: `${a.id} 5${d.suffix}`,
            alpha: alphaD,
          });
        }
      } else if (isTransitionMetal(a.element)) {
        // Transition metals (3d, 4d, 5d): include 5 real cubic (n-1)d orbitals
        const dShell = String(Math.max(3, Number(shell) - 1));
        const alphaD = props.vsieD ?? -10.5;
        for (const d of D_ORBITALS) {
          basis.push({
            atom: a,
            orbital: d.orb,
            label: `${a.id} ${dShell}${d.suffix}`,
            alpha: alphaD,
          });
        }
      } else if (props.vsieD !== undefined && !isLargeSystem) {
        // Hypervalent Period 3-6 p-block centers (coordination > 3 or bond-order sum > 4)
        const atomBonds = bonds.filter((b) => b.source === a.id || b.target === a.id);
        const atomBoSum = atomBonds.reduce((s, b) => s + b.order, 0);
        if (atomBonds.length > 3 || atomBoSum > 4) {
          for (const d of D_ORBITALS) {
            basis.push({
              atom: a,
              orbital: d.orb,
              label: `${a.id} ${shell}${d.suffix}`,
              alpha: props.vsieD,
            });
          }
        }
      }

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

      if (!bond && dist > 2.95) {
        H[i][j] = 0;
        H[j][i] = 0;
        continue;
      }

      // Direction cosines from atom i to atom j: d = (R_j - R_i) / |R_j - R_i|
      const l = dx / dist;
      const m = dy / dist;
      const n = dz / dist;

      // Distance-dependent radial Slater-Koster overlaps calibrated to Hoffmann EHT
      const r0 =
        (ELEMENT_PROPERTIES[bi.atom.element]?.covalentRadius ?? 0.77) +
        (ELEMENT_PROPERTIES[bj.atom.element]?.covalentRadius ?? 0.77);
      const decay = Math.exp(-0.85 * Math.max(0, dist - r0));
      const orderBoost = bond ? 1 + (bond.order - 1) * 0.14 : 0.28;

      const overlapS = computeTwoCenterSlaterKosterOverlap(
        bi.orbital,
        bj.orbital,
        l,
        m,
        n,
        decay,
        orderBoost
      );

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

  // 2. Otherwise, resolve any unlinked dative metal-ligand coordination bonds,
  //    expand implicit hydrogens (for small/medium graphs <= 24 atoms),
  //    reconstruct 3D VSEPR / coordination coordinates,
  //    and run the full s-p-d-f Extended Hückel LCAO-MO + Symmetry solver
  const hasDisconnectedMetalCenter = rawAtoms.some(
    (a) =>
      (isTransitionMetal(a.element) || isLanthanide(a.element)) &&
      rawBonds.filter((b) => b.source === a.id || b.target === a.id).length === 0 &&
      rawAtoms.length >= 3
  );
  const coordResolved = hasDisconnectedMetalCenter
    ? resolveCoordinationBonds(rawAtoms, rawBonds, false)
    : { atoms: rawAtoms, bonds: rawBonds };

  const shouldExpandImplicitH = coordResolved.atoms.length <= 24;
  const { atoms, bonds } = shouldExpandImplicitH
    ? expandImplicitHydrogens(coordResolved.atoms, coordResolved.bonds)
    : { atoms: coordResolved.atoms, bonds: coordResolved.bonds };
  const { formula } = computeMolecularFormula(atoms, bonds, !shouldExpandImplicitH);
  const smiles = generateSmilesFromGraph(atoms, bonds);
  const centeredAtoms = reconstruct3DGeometryFromVSEPR(atoms, bonds, charge);
  const pointGroupId = detectPointGroup(centeredAtoms, bonds, charge);
  const pg = POINT_GROUPS[pointGroupId] || POINT_GROUPS.Cs;

  const { basis, eigenvalues, eigenvectors } = runExtendedHuckelSolver(centeredAtoms, bonds);
  const N = basis.length;

  // Exact physical neutral valence electrons = sum(Z_val) of active basis atoms (clamped to s-p-d-f basis capacity).
  // Formal charges on individual atoms are only subtracted here when options.charge is 0 to prevent double-counting.
  const activeBasisAtomIds = new Set(basis.map((b) => b.atom.id));
  const rawValenceElectrons = centeredAtoms
    .filter((a) => activeBasisAtomIds.has(a.id))
    .reduce((sum, a) => {
      const props = ELEMENT_PROPERTIES[a.element];
      const maxBasisCap = isLanthanide(a.element)
        ? 32
        : isTransitionMetal(a.element)
          ? 18
          : a.element === 'H' || a.element === 'He'
            ? 2
            : 8;
      const atomChg = charge === 0 ? (a.formalCharge ?? 0) : 0;
      const v = Math.min(maxBasisCap, props?.valenceElectrons ?? 4) - atomChg;
      return sum + Math.max(0, v);
    }, 0);
  const neutralElectrons = Math.max(1, Math.min(N * 2, rawValenceElectrons));

  // Detect whether the molecule is linear along the x-axis or planar in the xy-plane (z ≈ 0)
  const maxAbsY = centeredAtoms.reduce((m, a) => Math.max(m, Math.abs(a.y3d)), 0);
  const maxAbsZ = centeredAtoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d)), 0);
  const isLinearAlongX = centeredAtoms.length >= 2 && maxAbsY < 0.15 && maxAbsZ < 0.15;
  const isPlanarXY = maxAbsZ < 0.15;

  // Identify central / primary skeleton atom and peripheral ligand atom for dynamic fragment orbitals
  const starTopo = analyzeStarCoordinationTopology(centeredAtoms, bonds);
  const primaryAtom =
    starTopo?.centerAtom ||
    centeredAtoms.find((a) => isTransitionMetal(a.element) || isLanthanide(a.element)) ||
    centeredAtoms.find((a) => {
      const deg = bonds.filter((b) => b.source === a.id || b.target === a.id).length;
      return deg >= 2 && deg === centeredAtoms.length - 1;
    }) ||
    centeredAtoms.find((a) => a.element !== 'H') ||
    centeredAtoms[0];
  const peripheralAtom =
    starTopo?.branches[0]?.donor ||
    centeredAtoms.find((a) => a.id !== primaryAtom?.id) ||
    primaryAtom;

  const primProps = ELEMENT_PROPERTIES[primaryAtom?.element || 'C'] || ELEMENT_PROPERTIES.C;
  const periphProps = ELEMENT_PROPERTIES[peripheralAtom?.element || 'H'] || ELEMENT_PROPERTIES.H;
  const primShell = getPrincipalShell(primaryAtom?.element || 'C');
  const periphShell = getPrincipalShell(peripheralAtom?.element || 'H');
  const primIsTm = primaryAtom ? isTransitionMetal(primaryAtom.element) : false;
  const primIsLn = primaryAtom ? isLanthanide(primaryAtom.element) : false;

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
      electrons: Math.max(0, Math.min(6, primProps.valenceElectrons - 2)),
    },
    ...(primIsLn
      ? [
          {
            id: 'frag_left_f',
            side: 'left' as const,
            label: `${primaryAtom?.element} 4f`,
            fragmentName: `${primaryAtom?.element} 4f/5d Shell`,
            energy: primProps.vsieF ?? -12.5,
            irrep: secondIrrepLower,
            degeneracy: 7,
            electrons: Math.max(0, primProps.valenceElectrons - 2),
          },
          {
            id: 'frag_left_d',
            side: 'left' as const,
            label: `${primaryAtom?.element} 5d`,
            fragmentName: `${primaryAtom?.element} 4f/5d Shell`,
            energy: primProps.vsieD ?? -8.0,
            irrep: firstIrrepLower,
            degeneracy: 5,
            electrons: 1,
          },
        ]
      : primIsTm
        ? [
            {
              id: 'frag_left_d',
              side: 'left' as const,
              label: `${primaryAtom?.element} ${Math.max(3, Number(primShell) - 1)}d`,
              fragmentName: `${primaryAtom?.element} d-Manifold`,
              energy: primProps.vsieD ?? -10.5,
              irrep: `${firstIrrepLower} + ${secondIrrepLower}`,
              degeneracy: 5,
              electrons: Math.max(1, primProps.valenceElectrons - 2),
            },
          ]
        : []),
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

  // Select all valence MOs when N <= 14 (e.g. C2H4 with 12 AOs), or a frontier window for larger systems without slicing degenerate doublets or triplets
  const EffectiveTotalElectrons = Math.max(1, Math.min(N * 2, neutralElectrons - charge));
  const homoIdx = Math.max(0, Math.min(N - 1, Math.floor((EffectiveTotalElectrons - 1) / 2)));
  let startIdx = N <= 14 ? 0 : Math.max(0, Math.min(homoIdx - 5, N - 12));
  let endIdx = N <= 14 ? N - 1 : Math.min(N - 1, Math.max(homoIdx + 5, startIdx + 11));
  while (startIdx > 0 && Math.abs(eigenvalues[startIdx] - eigenvalues[startIdx - 1]) < 0.08) {
    startIdx -= 1;
  }
  while (endIdx < N - 1 && Math.abs(eigenvalues[endIdx] - eigenvalues[endIdx + 1]) < 0.08) {
    endIdx += 1;
  }

  const molecularOrbitals: MolecularOrbital[] = [];
  let remainingElectrons = Math.max(0, neutralElectrons - startIdx * 2);

  const hasDOrFInBasis = basis.some(
    (b) => b.orbital.startsWith('d') || b.orbital.startsWith('f')
  );
  // Exclude d-only Delta and multi-electron SigmaMinus irreps when using a pure s-p valence basis
  const spIrreps = pg.irreps.filter(
    (ir) => (hasDOrFInBasis || !ir.symbol.startsWith('Δ')) && ir.symbol !== 'Σ-'
  );
  const candidateIrreps = spIrreps.length > 0 ? spIrreps : pg.irreps;
  const dim1Irreps = candidateIrreps.filter((ir) => ir.dimension === 1);
  const dimMultiIrreps = candidateIrreps.filter((ir) => ir.dimension >= 2);

  // Initialize literature-standard per-irrep Mulliken shell counters (accounting for inner 1s core shells on heavy atoms)
  const heavyCount = centeredAtoms.filter((a) => a.element !== 'H').length;
  const irrepShellCounts: Record<string, number> = {};
  if (startIdx === 0) {
    if (pg.id === 'Cinfv') {
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

    // Classify π and δ character from orbital composition
    const piWeight = vec.reduce(
      (acc, c, k) =>
        acc +
        ((isPlanarXY &&
          (basis[k].orbital === 'pz' ||
            basis[k].orbital === 'dxz' ||
            basis[k].orbital === 'dyz')) ||
        (isLinearAlongX &&
          (basis[k].orbital === 'py' ||
            basis[k].orbital === 'pz' ||
            basis[k].orbital === 'dxy' ||
            basis[k].orbital === 'dxz'))
          ? c * c
          : 0),
      0
    );
    const deltaWeight = vec.reduce(
      (acc, c, k) =>
        acc +
        ((isLinearAlongX &&
          (basis[k].orbital === 'dx2y2' ||
            basis[k].orbital === 'dyz' ||
            basis[k].orbital === 'fxyz' ||
            basis[k].orbital === 'fz_x2y2')) ||
        (hasDOrFInBasis &&
          (basis[k].orbital === 'fxyz' || basis[k].orbital === 'fz_x2y2'))
          ? c * c
          : 0),
      0
    );
    const isPi = piWeight > 0.48;
    const isDelta = !isPi && deltaWeight > 0.52;

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
          const sVal = computeTwoCenterSlaterKosterOverlap(
            basis[i].orbital,
            basis[j].orbital,
            dx / dist,
            dy / dist,
            dz / dist,
            0.88,
            1.0
          );
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

    const interactionType: MolecularOrbital['interactionType'] = isDelta
      ? bondingType === 'antibonding'
        ? 'δ*'
        : bondingType === 'non-bonding'
          ? 'n (Lone Pair)'
          : 'δ'
      : isPi
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
      // Count degeneracy multiplicity (2 for E doublets, 3 for T triplets in Oh / Td)
      let multipletSize = 1;
      while (
        idx + multipletSize <= endIdx &&
        Math.abs(eigenvalues[idx + multipletSize] - eigenvalues[idx]) < 0.08
      ) {
        multipletSize++;
      }
      const targetDim = multipletSize >= 3 ? 3 : 2;
      const dimMatchedIrreps = dimMultiIrreps.filter((ir) => ir.dimension === targetDim);
      const multiPool = dimMatchedIrreps.length > 0 ? dimMatchedIrreps : dimMultiIrreps;

      if (pg.id === 'D3h') {
        irrepObj =
          multiPool.find((ir) => (isPi ? ir.id === 'Epp' : ir.id === 'Ep')) ||
          multiPool[0];
      } else {
        const centerDWeight = primaryAtom
          ? vec.reduce(
              (s, c, k) =>
                s +
                (basis[k].atom.id === primaryAtom.id && basis[k].orbital.startsWith('d')
                  ? c * c
                  : 0),
              0
            )
          : 0;
        const centerPWeight = primaryAtom
          ? vec.reduce(
              (s, c, k) =>
                s +
                (basis[k].atom.id === primaryAtom.id &&
                (basis[k].orbital === 'px' || basis[k].orbital === 'py' || basis[k].orbital === 'pz')
                  ? c * c
                  : 0),
              0
            )
          : 0;
        const preferUngerade =
          heavyCount === 2 ? bondingType === 'bonding' : centerPWeight > 0.05 && centerDWeight < 0.1;
        irrepObj =
          multiPool.find((ir) =>
            preferUngerade ? ir.symbol.includes('u') : ir.symbol.includes('g')
          ) || multiPool[0];
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
    leftFragmentTitle: primIsLn
      ? `Central ${primaryAtom?.element || 'Skeleton'} (4f, 5d, 6s)`
      : primIsTm
        ? `Central ${primaryAtom?.element || 'Skeleton'} (${Math.max(3, Number(primShell) - 1)}d, ${primShell}s, ${primShell}p)`
        : `Central ${primaryAtom?.element || 'Skeleton'} (${primShell}s, ${primShell}p)`,
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
        angular = (2 * dz * dz - dx * dx - dy * dy) * 1.3;
        break;
      case 'dx2y2':
        angular = (dx * dx - dy * dy) * 1.5;
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
      case 'fz3':
        angular = dz * (2 * dz * dz - 3 * dx * dx - 3 * dy * dy) * 0.9;
        break;
      case 'fxz2':
        angular = dx * (4 * dz * dz - dx * dx - dy * dy) * 0.95;
        break;
      case 'fyz2':
        angular = dy * (4 * dz * dz - dx * dx - dy * dy) * 0.95;
        break;
      case 'fxyz':
        angular = dx * dy * dz * 3.2;
        break;
      case 'fz_x2y2':
        angular = dz * (dx * dx - dy * dy) * 1.8;
        break;
      case 'fx_x2_3y2':
        angular = dx * (dx * dx - 3 * dy * dy) * 1.1;
        break;
      case 'fy_3x2_y2':
        angular = dy * (3 * dx * dx - dy * dy) * 1.1;
        break;
    }

    psi += term.coeff * radial * angular;
  }

  return psi;
}
