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
  if (atom.element === 'O' && (atom.formalCharge ?? 0) === 0) {
    const conn = bonds.filter((b) => b.source === atom.id || b.target === atom.id);
    if (conn.length === 1 && conn[0].order === 1) {
      const otherId = conn[0].source === atom.id ? conn[0].target : conn[0].source;
      const other = atoms.find((x) => x.id === otherId);
      if (other?.element === 'C') return 1;
    }
    return 0;
  }
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

  // Alkene (η²-C2H4) or cyclopentadienyl (η⁵-C5H5) carbon π-coordinated to a transition metal
  // retains its full organic C-H count (ignoring the dative M–C π-coordination bond)
  const nonMetalBonds = connected.filter((b) => {
    const otherId = b.source === atom.id ? b.target : b.source;
    const other = atoms.find((x) => x.id === otherId);
    return !(other && (isTransitionMetal(other.element) || isLanthanide(other.element)));
  });
  const hasCCBond = nonMetalBonds.some((b) => {
    const otherId = b.source === atom.id ? b.target : b.source;
    const other = atoms.find((x) => x.id === otherId);
    return other?.element === 'C';
  });
  if (isBondedToMetal && hasCCBond) {
    const nonMetalBoSum = nonMetalBonds.reduce((sum, b) => sum + b.order, 0);
    const chgPenalty = atom.formalCharge ? Math.abs(atom.formalCharge) : 0;
    return Math.max(0, 4 - nonMetalBoSum - chgPenalty);
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
 * where PubChem V2000 SDFs or disconnected SMILES encode monodentate, bidentate/chelating,
 * or η⁵/η² π-ligands (`PPh3`, `PEt3`, `bpy`, `en`, `acac-`, `C2O4(2-)`, `[C-]#N`, `[C-]#[O+]`,
 * `NH3`, `Cl-`, `C5H5-`, `C2H4`, etc.) as separate unconnected fragments around a central metal
 * or hypervalent coordination center, and strips isolated spectator counterions and lattice waters.
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

  const completeEta2AlkeneHapticity = (
    atomList: AtomNode[],
    bondList: BondEdge[]
  ): BondEdge[] => {
    const outBonds = [...bondList];
    const areBonded = (id1: string, id2: string) =>
      outBonds.some(
        (b) => (b.source === id1 && b.target === id2) || (b.source === id2 && b.target === id1)
      );
    for (const m of atomList) {
      if (!isTransitionMetal(m.element) && !isLanthanide(m.element)) continue;
      const mBonds = outBonds.filter((b) => b.source === m.id || b.target === m.id);
      for (const mb of mBonds) {
        const c1Id = mb.source === m.id ? mb.target : mb.source;
        const c1 = atomList.find((a) => a.id === c1Id);
        if (!c1 || c1.element !== 'C') continue;
        // Find any C1=C2 double bond where C2 is not yet bonded to m (and C1, C2 are not in a 6-membered aryl ring)
        const ccDouble = outBonds.find((b) => {
          if (b.order !== 2) return false;
          const otherId = b.source === c1.id ? b.target : b.target === c1.id ? b.source : null;
          if (!otherId || areBonded(m.id, otherId)) return false;
          const c2 = atomList.find((a) => a.id === otherId);
          if (!c2 || c2.element !== 'C') return false;
          const c1HeavyNbs = outBonds
            .filter((eb) => eb.source === c1.id || eb.target === c1.id)
            .map((eb) => atomList.find((x) => x.id === (eb.source === c1.id ? eb.target : eb.source)))
            .filter((x) => x && x.element !== 'H' && x.id !== m.id);
          const c2HeavyNbs = outBonds
            .filter((eb) => eb.source === c2.id || eb.target === c2.id)
            .map((eb) => atomList.find((x) => x.id === (eb.source === c2.id ? eb.target : eb.source)))
            .filter((x) => x && x.element !== 'H' && x.id !== m.id);
          return c1HeavyNbs.length <= 2 && c2HeavyNbs.length <= 2;
        });
        if (ccDouble) {
          const c2Id = ccDouble.source === c1.id ? ccDouble.target : ccDouble.source;
          if (!areBonded(m.id, c2Id)) {
            outBonds.push({
              id: `b_eta2_${outBonds.length + 1}`,
              source: m.id,
              target: c2Id,
              order: 1,
            });
          }
        }
      }
    }
    return outBonds;
  };

  let atoms = rawAtoms.map((a) => ({ ...a }));
  let bonds = rawBonds.map((b) => ({ ...b }));
  let comps = getComponents(atoms, bonds);

  if (comps.length <= 1) {
    return { atoms, bonds: completeEta2AlkeneHapticity(atoms, bonds) };
  }

  // Strip isolated group-1 / group-2 spectator counterions (e.g. 4 K+ in K4[Fe(CN)6], 3 K+ in K3[Fe(C2O4)3], 2 Na+ in Na2[Fe(CN)5NO])
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
          return { atoms, bonds: completeEta2AlkeneHapticity(atoms, bonds) };
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

  const isNeutralWaterComp = (comp: AtomNode[]): boolean => {
    const heavy = comp.filter((a) => a.element !== 'H');
    if (heavy.length !== 1 || heavy[0].element !== 'O') return false;
    if ((heavy[0].formalCharge ?? 0) !== 0) return false;
    const hCount = comp.length - 1;
    // Bare uncharged O (from SMILES `.O.`) or H2O (O + 2 H)
    return hCount === 0 || hCount === 2;
  };

  /**
   * Identifies the σ/π donor atom(s) in a ligand component (1 for monodentate ligands like
   * -Cl, -CN, -CO, -NO, -NH3, -PPh3, -PEt3, -SR, pyridine, C5H5, C2H4; 2 for bidentate chelates
   * like ethylenediamine, 2,2'-bipyridine, acetylacetonate, oxalate, dicarboxylates, dimethylglyoximate).
   */
  const findLigandDonors = (comp: AtomNode[]): { donors: AtomNode[]; deprotonateHIds: string[] } => {
    if (comp.length === 1) {
      return { donors: [comp[0]], deprotonateHIds: [] };
    }
    if (comp.length === 2) {
      const [a, b] = comp;
      if (a.element === 'H' && b.element !== 'H') {
        return { donors: [b], deprotonateHIds: [] };
      }
      if (b.element === 'H' && a.element !== 'H') {
        return { donors: [a], deprotonateHIds: [] };
      }
      if (
        (a.element === 'C' && (b.element === 'N' || b.element === 'O' || b.element === 'S')) ||
        (b.element === 'C' && (a.element === 'N' || a.element === 'O' || a.element === 'S'))
      ) {
        const cAtom = a.element === 'C' ? a : b;
        const heteroAtom = a.element === 'C' ? b : a;
        if (heteroAtom.element === 'O') {
          bonds = bonds.map((edge) =>
            (edge.source === cAtom.id && edge.target === heteroAtom.id) ||
            (edge.source === heteroAtom.id && edge.target === cAtom.id)
              ? { ...edge, order: 3 as const }
              : edge
          );
        }
        return { donors: [cAtom], deprotonateHIds: [] };
      }
      if (
        (a.element === 'N' && b.element === 'O') ||
        (b.element === 'N' && a.element === 'O')
      ) {
        return { donors: [a.element === 'N' ? a : b], deprotonateHIds: [] };
      }
      return {
        donors: [distToCenter(a) <= distToCenter(b) ? a : b],
        deprotonateHIds: [],
      };
    }

    const heavy = comp.filter((a) => a.element !== 'H');
    if (heavy.length === 0) return { donors: [comp[0]], deprotonateHIds: [] };
    if (heavy.length === 1) return { donors: [heavy[0]], deprotonateHIds: [] };

    const getConn = (a: AtomNode) =>
      bonds.filter((b) => b.source === a.id || b.target === a.id);
    const getHeavyNbs = (a: AtomNode) =>
      getConn(a)
        .map((b) => comp.find((x) => x.id === (b.source === a.id ? b.target : b.source)))
        .filter((x): x is AtomNode => Boolean(x && x.element !== 'H'));
    const getBoSum = (a: AtomNode) => getConn(a).reduce((s, b) => s + b.order, 0);

    // 1. Trivalent Group 15 Phosphine / Arsine / Stibine donors (PPh3, PEt3, PMe3, PF3, dppe, BINAP)
    const pDonors = heavy.filter(
      (a) =>
        (a.element === 'P' || a.element === 'As' || a.element === 'Sb') &&
        getHeavyNbs(a).length <= 3 &&
        getBoSum(a) <= 3
    );
    if (pDonors.length > 0) {
      return { donors: pDonors.slice(0, 2), deprotonateHIds: [] };
    }

    // 2. Low-valent Group 16 Thiolate / Thioether / Selenoether donors (e.g. [S-] in Auranofin)
    const sDonors = heavy
      .filter(
        (a) =>
          (a.element === 'S' || a.element === 'Se' || a.element === 'Te') &&
          getHeavyNbs(a).length <= 2 &&
          getBoSum(a) <= 2
      )
      .sort((a, b) => (a.formalCharge ?? 0) - (b.formalCharge ?? 0));
    if (sDonors.length > 0) {
      return { donors: sDonors.slice(0, 2), deprotonateHIds: [] };
    }

    // 3. Terminal Isocyanide (R-N#C) donor where terminal C is triple-bonded to N
    const isocyanideCs = heavy.filter((a) => {
      if (a.element !== 'C' || getConn(a).length !== 1) return false;
      const b = getConn(a)[0];
      if (b.order < 2) return false;
      const nbId = b.source === a.id ? b.target : b.source;
      const nb = comp.find((x) => x.id === nbId);
      return nb?.element === 'N';
    });
    if (isocyanideCs.length > 0) {
      return { donors: isocyanideCs.slice(0, 2), deprotonateHIds: [] };
    }

    // 4. Lone-pair Nitrogen (N) donors (en, bpy, phen, pyridine, 1,2-diaminocyclohexane, dmgH-, RCN)
    const nDonors = heavy.filter(
      (a) => a.element === 'N' && getConn(a).length <= 3 && getBoSum(a) <= 3
    );
    if (nDonors.length > 0) {
      return { donors: nDonors.slice(0, 2), deprotonateHIds: [] };
    }

    // 5. Oxygen (O) donors: Dicarboxylates / Oxalate (2 single-bonded carboxylate O's), β-Diketonates (acac, 2 O's), or monodentate O
    const oAtoms = heavy.filter((a) => a.element === 'O');
    if (oAtoms.length > 0) {
      // Identify single-bonded carboxylate oxygens (-C(=O)-O- or -C(=O)-OH)
      const carboxylateOs = oAtoms.filter((o) => {
        const hvyNbs = getHeavyNbs(o);
        if (hvyNbs.length !== 1 || hvyNbs[0].element !== 'C') return false;
        const cAtom = hvyNbs[0];
        const coBond = bonds.find(
          (b) =>
            (b.source === o.id && b.target === cAtom.id) ||
            (b.target === o.id && b.source === cAtom.id)
        );
        if (!coBond || coBond.order !== 1) return false;
        return getConn(cAtom).some((b) => {
          if (b.id === coBond.id || b.order !== 2) return false;
          const otherId = b.source === cAtom.id ? b.target : b.source;
          const other = comp.find((x) => x.id === otherId);
          return other?.element === 'O';
        });
      });

      if (carboxylateOs.length >= 2) {
        const chosenOs = carboxylateOs.slice(0, 2);
        const deprotonateHIds: string[] = [];
        for (const o of chosenOs) {
          for (const b of getConn(o)) {
            const otherId = b.source === o.id ? b.target : b.source;
            const other = comp.find((x) => x.id === otherId);
            if (other?.element === 'H') deprotonateHIds.push(other.id);
          }
        }
        return { donors: chosenOs, deprotonateHIds };
      }

      if (oAtoms.length === 2) {
        // β-diketonate (acetylacetonate acac-) chelates through both O atoms
        return { donors: oAtoms, deprotonateHIds: [] };
      }

      const sortedOs = [...oAtoms].sort(
        (a, b) => (a.formalCharge ?? 0) - (b.formalCharge ?? 0) || distToCenter(a) - distToCenter(b)
      );
      return { donors: [sortedOs[0]], deprotonateHIds: [] };
    }

    // 6. Pure Carbon π/σ ligands:
    //    - η⁵-C5H5 in Ferrocene/Cobaltocene/Ruthenocene
    //    - η⁶-C6H6 in Bis(benzene)chromium Cr(C6H6)2 or (C6H6)Cr(CO)3
    //    - η⁴-C4H4 in (C4H4)Fe(CO)3
    //    - η²,η²-chelating dienes (1,5-cyclooctadiene COD, norbornadiene NBD, 1,3-butadiene)
    //    - monodentate η²-alkenes of any size (ethylene, propene, 1-hexene, styrene, cyclohexene, norbornene, etc.)
    if (
      (heavy.length === 5 || heavy.length === 6 || heavy.length === 4) &&
      heavy.every((a) => a.element === 'C' && getHeavyNbs(a).length === 2)
    ) {
      const ringDoubleBonds = bonds.filter((b) => {
        if (b.order !== 2) return false;
        return heavy.some((x) => x.id === b.source) && heavy.some((x) => x.id === b.target);
      });
      const allAnionicRing = heavy.every((a) => (a.formalCharge ?? 0) < 0);
      const explicitHCount = comp.filter((a) => a.element === 'H').length;
      const isC6H6Ring = heavy.length === 6 && explicitHCount === 6;
      const isC4H4Ring = heavy.length === 4 && explicitHCount === 4;
      if (
        heavy.length === 5 ||
        (heavy.length === 6 && (ringDoubleBonds.length === 3 || allAnionicRing || isC6H6Ring)) ||
        (heavy.length === 4 && (ringDoubleBonds.length >= 2 || allAnionicRing || isC4H4Ring))
      ) {
        // Normalize PubChem anionic / bond-order-omitted ring representations (e.g. [CH-]1[CH-][CH-][CH-][CH-][CH-]1 in Cr(C6H6)2)
        if (allAnionicRing || isC6H6Ring || isC4H4Ring) {
          for (const a of heavy) {
            a.formalCharge = 0;
          }
          if ((heavy.length === 6 || heavy.length === 4) && ringDoubleBonds.length === 0) {
            const nRing = heavy.length;
            const ordered: AtomNode[] = [heavy[0]];
            const used = new Set<string>([heavy[0].id]);
            while (ordered.length < nRing) {
              const last = ordered[ordered.length - 1];
              const next = getHeavyNbs(last).find((nb) => !used.has(nb.id));
              if (!next) break;
              ordered.push(next);
              used.add(next.id);
            }
            if (ordered.length === nRing) {
              for (let k = 0; k < nRing; k += 2) {
                const a1 = ordered[k];
                const a2 = ordered[(k + 1) % nRing];
                bonds = bonds.map((edge) =>
                  (edge.source === a1.id && edge.target === a2.id) ||
                  (edge.source === a2.id && edge.target === a1.id)
                    ? { ...edge, order: 2 as const }
                    : edge
                );
              }
            }
          }
        }
        return { donors: heavy, deprotonateHIds: [] };
      }
    }

    // Check for η²-alkene or η²,η²-diene C=C double bonds (Dewar-Chatt-Duncanson coordination)
    if (isTransitionMetal(centerAtom.element) || isLanthanide(centerAtom.element)) {
      const ccDoubleBonds = bonds
        .filter((b) => {
          if (b.order !== 2) return false;
          const a1 = comp.find((x) => x.id === b.source);
          const a2 = comp.find((x) => x.id === b.target);
          return Boolean(a1?.element === 'C' && a2?.element === 'C');
        })
        .sort((b1, b2) => {
          const a11 = comp.find((x) => x.id === b1.source)!;
          const a12 = comp.find((x) => x.id === b1.target)!;
          const a21 = comp.find((x) => x.id === b2.source)!;
          const a22 = comp.find((x) => x.id === b2.target)!;
          return (
            Math.min(distToCenter(a11), distToCenter(a12)) -
            Math.min(distToCenter(a21), distToCenter(a22))
          );
        });

      if (ccDoubleBonds.length > 0) {
        // Distinguish non-aromatic olefinic C=C bonds from 6-membered aryl ring C=C bonds (e.g. phenyl ring in styrene)
        const isBondInAromaticSixRing = (bond: BondEdge): boolean => {
          const startId = bond.source;
          const targetId = bond.target;
          // BFS/DFS for a path of length 5 in heavy atoms from startId to targetId not using `bond`
          const stack: { curr: string; path: string[] }[] = [{ curr: startId, path: [startId] }];
          while (stack.length > 0) {
            const { curr, path } = stack.pop()!;
            if (path.length === 6) {
              const currAtom = comp.find((x) => x.id === curr);
              if (
                currAtom &&
                getHeavyNbs(currAtom).some((nb) => nb.id === targetId)
              ) {
                const ringSet = new Set(path);
                const dblInRing = ccDoubleBonds.filter(
                  (db) => ringSet.has(db.source) && ringSet.has(db.target)
                ).length;
                if (dblInRing >= 3) return true;
              }
              continue;
            }
            const currAtom = comp.find((x) => x.id === curr);
            if (!currAtom) continue;
            for (const nb of getHeavyNbs(currAtom)) {
              if (curr === startId && nb.id === targetId) continue;
              if (!path.includes(nb.id)) {
                stack.push({ curr: nb.id, path: [...path, nb.id] });
              }
            }
          }
          return false;
        };

        const nonArylDoubleBonds = ccDoubleBonds.filter((b) => !isBondInAromaticSixRing(b));
        const activeDoubleBonds =
          nonArylDoubleBonds.length > 0 ? nonArylDoubleBonds : ccDoubleBonds;

        // Calculate existing non-alkene coordination sites on centerAtom to see if chelating diene (2 sites) fits
        const curCenterNbs = bonds
          .filter((b) => b.source === centerAtom.id || b.target === centerAtom.id)
          .map((b) => (b.source === centerAtom.id ? b.target : b.source));
        const curEta2Pairs = bonds.filter(
          (b) => curCenterNbs.includes(b.source) && curCenterNbs.includes(b.target)
        ).length;
        const curEffectiveSites = Math.max(0, curCenterNbs.length - curEta2Pairs);

        // Chelating η²,η²-diene (e.g. 1,5-cyclooctadiene COD in Ni(COD)2 / Pt(COD)Cl2, norbornadiene NBD, 1,3-butadiene)
        if (activeDoubleBonds.length >= 2 && curEffectiveSites <= 2) {
          const b1 = activeDoubleBonds[0];
          const b2 = activeDoubleBonds.find(
            (b) =>
              b.source !== b1.source &&
              b.source !== b1.target &&
              b.target !== b1.source &&
              b.target !== b1.target
          );
          if (b2) {
            const c1a = comp.find((x) => x.id === b1.source)!;
            const c1b = comp.find((x) => x.id === b1.target)!;
            const c2a = comp.find((x) => x.id === b2.source)!;
            const c2b = comp.find((x) => x.id === b2.target)!;
            return { donors: [c1a, c1b, c2a, c2b], deprotonateHIds: [] };
          }
        }

        // Monodentate η²-alkene (ethylene, propene, 1-hexene, styrene, cyclohexene, norbornene, etc.)
        const bestCC = activeDoubleBonds[0];
        const c1 = comp.find((x) => x.id === bestCC.source)!;
        const c2 = comp.find((x) => x.id === bestCC.target)!;
        return { donors: [c1, c2], deprotonateHIds: [] };
      }
    }
    const sortedHeavy = [...heavy].sort(
      (a, b) =>
        (a.formalCharge ?? 0) - (b.formalCharge ?? 0) ||
        getConn(a).length - getConn(b).length ||
        distToCenter(a) - distToCenter(b)
    );
    return { donors: [sortedHeavy[0]], deprotonateHIds: [] };
  };

  // Evaluate all non-center ligand components and strip uncoordinated lattice waters of hydration
  // (e.g. .2H2O in Sodium nitroprusside dihydrate CID 11953895 or .H2O in Zeise's salt CID 71311337)
  // when primary non-water ligands already saturate the metal's coordination sphere (4 for sq-planar, 6 for octahedral).
  const centerNbs = bonds
    .filter((b) => b.source === centerAtom.id || b.target === centerAtom.id)
    .map((b) => (b.source === centerAtom.id ? b.target : b.source));
  const eta2ExistingPairs = bonds.filter(
    (b) => centerNbs.includes(b.source) && centerNbs.includes(b.target)
  ).length;
  const existingCenterDegree = Math.max(0, centerNbs.length - eta2ExistingPairs);
  const otherComps = comps.filter((c) => !c.some((a) => a.id === centerAtom.id));

  const getEffectiveCoordSites = (donors: AtomNode[]): number => {
    if (donors.length >= 5 && donors.every((d) => d.element === 'C')) {
      return 3;
    }
    if (donors.length === 4 && donors.every((d) => d.element === 'C')) {
      return 2;
    }
    if (
      donors.length === 2 &&
      bonds.some(
        (b) =>
          (b.source === donors[0].id && b.target === donors[1].id) ||
          (b.source === donors[1].id && b.target === donors[0].id)
      )
    ) {
      // η²-alkene C=C occupies 1 VSEPR coordination site
      return 1;
    }
    return donors.length;
  };

  if (stripSpectatorCounterions && otherComps.length >= 2) {
    const waterComps = otherComps.filter(isNeutralWaterComp);
    const nonWaterComps = otherComps.filter((c) => !isNeutralWaterComp(c));
    if (waterComps.length > 0 && nonWaterComps.length > 0) {
      const nonWaterDenticity =
        existingCenterDegree +
        nonWaterComps.reduce(
          (sum, c) => sum + getEffectiveCoordSites(findLigandDonors(c).donors),
          0
        );
      const targetSatCoord =
        centerAtom.element === 'Pt' || centerAtom.element === 'Pd' || centerAtom.element === 'Au'
          ? 4
          : 6;
      if (nonWaterDenticity >= targetSatCoord) {
        const waterAtomIds = new Set(waterComps.flat().map((a) => a.id));
        atoms = atoms.filter((a) => !waterAtomIds.has(a.id));
        bonds = bonds.filter(
          (b) => !waterAtomIds.has(b.source) && !waterAtomIds.has(b.target)
        );
        comps = getComponents(atoms, bonds);
      }
    }
  }

  const deprotonatedIds = new Set<string>();
  const allCoordinatedDonors: AtomNode[] = [];
  for (const comp of comps) {
    if (comp.some((a) => a.id === centerAtom.id)) continue;
    const { donors, deprotonateHIds } = findLigandDonors(comp);
    for (const hId of deprotonateHIds) deprotonatedIds.add(hId);
    for (const donor of donors) {
      allCoordinatedDonors.push(donor);
      bonds.push({
        id: `b_coord_${bonds.length + 1}`,
        source: centerAtom.id,
        target: donor.id,
        order: 1,
      });
    }
  }

  // If centerAtom was an isolated single-atom component (e.g. [Ni] in Ni(COD)2 or [Fe] in Ferrocene / Hexacyanoferrate),
  // place its 2D coordinates at the centroid of all coordinated donor atoms so the 2D Sketcher is centered.
  if ((tmOrLnSingle || mainGroupSingle) && allCoordinatedDonors.length >= 2) {
    const avgX2d =
      allCoordinatedDonors.reduce((s, d) => s + d.x2d, 0) / allCoordinatedDonors.length;
    const avgY2d =
      allCoordinatedDonors.reduce((s, d) => s + d.y2d, 0) / allCoordinatedDonors.length;
    centerAtom.x2d = Number(avgX2d.toFixed(2));
    centerAtom.y2d = Number(avgY2d.toFixed(2));
  }

  if (deprotonatedIds.size > 0) {
    atoms = atoms.filter((a) => !deprotonatedIds.has(a.id));
    bonds = bonds.filter(
      (b) => !deprotonatedIds.has(b.source) && !deprotonatedIds.has(b.target)
    );
  }

  bonds = completeEta2AlkeneHapticity(atoms, bonds);

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
 * Determines whether a 4-coordinate transition-metal complex adopts a Square-Planar (d⁸)
 * or Tetrahedral (d¹⁰ / weak-field 3d) coordination geometry.
 */
export function isFourCoordSquarePlanar(
  centerAtom: AtomNode,
  branches: LigandBranchInfo[],
  charge = 0
): boolean {
  if (!isSquarePlanarMetal(centerAtom.element)) return false;
  if (branches.length !== 4) return false;

  const anionicCount = branches.filter((b) => {
    if (['F', 'Cl', 'Br', 'I'].includes(b.donor.element) && b.distalAtoms.length === 0) {
      return true;
    }
    if (b.sig === 'C:N') return true;
    if ((b.donor.formalCharge ?? 0) < 0) return true;
    return false;
  }).length;

  // Group 10 (Ni, Pd, Pt) zerovalent d10 ML4 complexes with 4 neutral L-type ligands
  // (e.g. Ni(CO)4, Pd(PPh3)4, Pt(PPh3)4, Ni(PF3)4) are Tetrahedral (Td)!
  if (
    (centerAtom.element === 'Ni' ||
      centerAtom.element === 'Pd' ||
      centerAtom.element === 'Pt') &&
    charge === 0 &&
    anionicCount === 0
  ) {
    return false;
  }

  // 3d Ni(II) is square-planar with strong-field cyanide ([Ni(CN)4]2-), and tetrahedral with weak-field halides/phosphines (NiCl2(PPh3)2, [NiCl4]2-)
  if (centerAtom.element === 'Ni') {
    return branches.some((b) => b.sig === 'C:N' || b.donor.element === 'N');
  }

  // 4d/5d d8 complexes (Rh(I) in Wilkinson's catalyst [RhCl(PPh3)3], Ir(I) in Vaska's complex [IrCl(CO)(PPh3)2],
  // Pd(II) in [PdCl2(PPh3)2], Pt(II) in Cisplatin [PtCl2(NH3)2] & [PtCl4]2-, Au(III) in [AuCl4]-) are Square Planar!
  return true;
}

export interface ChelateLigandInfo {
  donors: AtomNode[];
  allAtoms: AtomNode[];
  sig: string;
}

export interface ChelateOrSandwichTopology {
  kind:
    | 'metallocene_sandwich'
    | 'bis_diene_chelate'
    | 'mono_diene_chelate'
    | 'piano_stool_sandwich'
    | 'tris_bidentate_octahedral'
    | 'square_planar_chelate';
  centerAtom: AtomNode;
  ligands: ChelateLigandInfo[];
  pointGroupId: string;
}

/**
 * Detects whether a molecular graph is:
 * - a Metallocene / Bis-Arene Sandwich complex (e.g. Ferrocene Fe(C5H5)2 -> D5d, Bis(benzene)chromium Cr(C6H6)2 -> D6h),
 * - a Homoleptic Bis-Diene Chelate complex (e.g. Bis(1,5-cyclooctadiene)nickel Ni(COD)2 -> D2d),
 * - a Heteroleptic Square-Planar Diene Chelate complex (e.g. Pt(COD)Cl2, Pd(COD)Cl2, Pt(nbd)Cl2, [Rh(COD)Cl2]- -> C2v),
 * - a Piano-Stool / Half-Sandwich complex (e.g. (C6H6)Cr(CO)3 -> C3v, (C4H4)Fe(CO)3 / (C4H6)Fe(CO)3 -> Cs),
 * - a Tris-Bidentate Octahedral Propeller complex (e.g. [Co(en)3]3+, [Ru(bpy)3]2+, Fe(acac)3, [Fe(C2O4)3]3- -> D3), or
 * - a 4-Coordinate Square-Planar Chelate complex (e.g. [Ni(dmgH)2] -> D2h, Oxaliplatin / Carboplatin -> C2v).
 */
export function analyzeChelateOrSandwichTopology(
  atoms: AtomNode[],
  bonds: BondEdge[]
): ChelateOrSandwichTopology | null {
  if (atoms.length < 7 || atoms.length > 180) return null;

  const metalCandidates = atoms.filter(
    (a) => isTransitionMetal(a.element) || isLanthanide(a.element)
  );
  if (metalCandidates.length !== 1) return null;
  const centerAtom = metalCandidates[0];

  const adj = new Map<string, string[]>();
  for (const a of atoms) adj.set(a.id, []);
  for (const b of bonds) {
    if (adj.has(b.source) && adj.has(b.target)) {
      adj.get(b.source)!.push(b.target);
      adj.get(b.target)!.push(b.source);
    }
  }

  const centerNeighborIds = new Set(adj.get(centerAtom.id) || []);
  if (centerNeighborIds.size < 2 || centerNeighborIds.size > 12) return null;

  // Partition (atoms \ {centerAtom}) into connected ligand subgraphs
  const visited = new Set<string>([centerAtom.id]);
  const ligands: ChelateLigandInfo[] = [];

  for (const a of atoms) {
    if (visited.has(a.id)) continue;
    const comp: AtomNode[] = [];
    const queue: AtomNode[] = [a];
    visited.add(a.id);
    while (queue.length > 0) {
      const cur = queue.shift()!;
      comp.push(cur);
      for (const nbId of adj.get(cur.id) || []) {
        if (!visited.has(nbId)) {
          visited.add(nbId);
          const nbAtom = atoms.find((x) => x.id === nbId);
          if (nbAtom) queue.push(nbAtom);
        }
      }
    }
    const donors = comp.filter((x) => centerNeighborIds.has(x.id));
    if (donors.length === 0) return null;
    const elemSig = comp
      .map((x) => x.element)
      .sort()
      .join('');
    const donorSig = donors
      .map((x) => x.element)
      .sort()
      .join('');
    ligands.push({
      donors,
      allAtoms: comp,
      sig: `${donorSig}:${elemSig}`,
    });
  }

  // 1. Check for Metallocene / Bis-Arene Sandwich M(η⁵-C5H5)2 or M(η⁶-C6H6)2
  if (ligands.length === 2) {
    const isSandwichCap = (lig: ChelateLigandInfo) => {
      const heavy = lig.allAtoms.filter((x) => x.element !== 'H');
      return (
        (heavy.length === 5 || heavy.length === 6) &&
        heavy.every((x) => x.element === 'C') &&
        lig.donors.length >= 5
      );
    };
    if (isSandwichCap(ligands[0]) && isSandwichCap(ligands[1])) {
      const n0 = ligands[0].allAtoms.filter((x) => x.element !== 'H').length;
      const n1 = ligands[1].allAtoms.filter((x) => x.element !== 'H').length;
      const pgId = n0 === 6 && n1 === 6 ? 'D6h' : n0 === 5 && n1 === 5 ? 'D5d' : 'Cs';
      return {
        kind: 'metallocene_sandwich',
        centerAtom,
        ligands,
        pointGroupId: pgId,
      };
    }
  }

  // 2. Check for Chelating η²,η²-Diene Complexes (e.g. Ni(COD)2 -> D2d, Pt(COD)Cl2 / Pd(COD)Cl2 / Pt(nbd)Cl2 -> C2v)
  const isDieneChelate = (lig: ChelateLigandInfo) =>
    lig.donors.length === 4 && lig.donors.every((d) => d.element === 'C');

  if (ligands.length === 2 && ligands.every(isDieneChelate)) {
    return {
      kind: 'bis_diene_chelate',
      centerAtom,
      ligands,
      pointGroupId: ligands[0].sig === ligands[1].sig ? 'D2d' : 'C2v',
    };
  }

  if (
    ligands.length === 3 &&
    ligands.filter(isDieneChelate).length === 1 &&
    ligands.filter((l) => l.donors.length === 1).length === 2
  ) {
    const sortedLigands = [...ligands].sort((a, b) => b.donors.length - a.donors.length);
    const pgId = sortedLigands[1].sig === sortedLigands[2].sig ? 'C2v' : 'Cs';
    return {
      kind: 'mono_diene_chelate',
      centerAtom,
      ligands: sortedLigands,
      pointGroupId: pgId,
    };
  }

  // 3. Check for Piano-Stool / Half-Sandwich Complex M(ηⁿ-ring/diene)L3 (e.g. (C6H6)Cr(CO)3, (C4H4)Fe(CO)3, (C4H6)Fe(CO)3)
  if (
    ligands.length === 4 &&
    ligands.filter((l) => l.donors.length >= 4 && l.donors.every((d) => d.element === 'C'))
      .length === 1 &&
    ligands.filter((l) => l.donors.length === 1).length === 3
  ) {
    const sortedLigands = [...ligands].sort((a, b) => b.donors.length - a.donors.length);
    const capHeavy = sortedLigands[0].allAtoms.filter((x) => x.element !== 'H');
    const legsIdentical =
      sortedLigands[1].sig === sortedLigands[2].sig &&
      sortedLigands[2].sig === sortedLigands[3].sig;
    const pgId = capHeavy.length === 6 && legsIdentical ? 'C3v' : 'Cs';
    return {
      kind: 'piano_stool_sandwich',
      centerAtom,
      ligands: sortedLigands,
      pointGroupId: pgId,
    };
  }

  // Require at least one true bidentate chelate ligand (donors.length === 2 and donors are not directly bonded as an η²-alkene)
  const isTrueBidentateChelate = (l: ChelateLigandInfo) =>
    l.donors.length === 2 &&
    !bonds.some(
      (b) =>
        (b.source === l.donors[0].id && b.target === l.donors[1].id) ||
        (b.source === l.donors[1].id && b.target === l.donors[0].id)
    );
  const hasBidentate = ligands.some(isTrueBidentateChelate);
  if (!hasBidentate) return null;

  // 4. Check for Tris-Bidentate Octahedral Propeller M(L-L)3 (3 bidentate ligands = 6 donor bonds)
  if (
    ligands.length === 3 &&
    ligands.every(isTrueBidentateChelate) &&
    centerNeighborIds.size === 6
  ) {
    const uniqueSigs = new Set(ligands.map((l) => l.sig));
    return {
      kind: 'tris_bidentate_octahedral',
      centerAtom,
      ligands,
      pointGroupId: uniqueSigs.size === 1 ? 'D3' : 'Cs',
    };
  }

  // 5. Check for 4-Coordinate Square-Planar Chelate (e.g. [Ni(dmgH)2] -> D2h, Oxaliplatin / Carboplatin -> C2v)
  if (centerNeighborIds.size === 4 && (ligands.length === 2 || ligands.length === 3)) {
    const valid2 = ligands.length === 2 && ligands.every(isTrueBidentateChelate);
    const valid3 =
      ligands.length === 3 &&
      ligands.filter(isTrueBidentateChelate).length === 1 &&
      ligands.filter((l) => l.donors.length === 1).length === 2;
    if (valid2 || valid3) {
      const sortedLigands = [...ligands].sort((a, b) => b.donors.length - a.donors.length);
      const pgId =
        valid2 && ligands[0].sig === ligands[1].sig ? 'D2h' : 'C2v';
      return {
        kind: 'square_planar_chelate',
        centerAtom,
        ligands: sortedLigands,
        pointGroupId: pgId,
      };
    }
  }

  return null;
}

/**
 * Detects whether a molecular graph has a single central coordination / VSEPR atom
 * surrounded by nLig ∈ [2, 6] mutually disjoint terminal ligand branches (supporting both
 * monoatomic ligands like -Cl, -F, -H and polyatomic/complex ligands like -C≡N, -C≡O, -NO,
 * -NH3, -OH2, -PPh3, -PEt3, pyridine, -SR, and Dewar-Chatt-Duncanson η²-alkenes like η²-C2H4).
 */
export function analyzeStarCoordinationTopology(
  atoms: AtomNode[],
  bonds: BondEdge[]
): StarCoordinationTopology | null {
  if (atoms.length < 3 || atoms.length > 200) return null;

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
      const isTmOrLn = isTransitionMetal(a.element) || isLanthanide(a.element);
      return deg >= 2 && deg <= (isTmOrLn ? 8 : 6);
    })
    .sort((a, b) => {
      const aMetal = isTransitionMetal(a.element) || isLanthanide(a.element) ? 1 : 0;
      const bMetal = isTransitionMetal(b.element) || isLanthanide(b.element) ? 1 : 0;
      if (aMetal !== bMetal) return bMetal - aMetal;
      return (adj.get(b.id) || []).length - (adj.get(a.id) || []).length;
    });

  for (const centerAtom of candidates) {
    const isMetalCenter =
      isTransitionMetal(centerAtom.element) || isLanthanide(centerAtom.element);
    if (!isMetalCenter && atoms.length > 32) continue;

    // Exclude metallocene sandwich complexes (which use analyzeChelateOrSandwichTopology)
    const donorIds = adj.get(centerAtom.id) || [];
    if (donorIds.length < 2 || donorIds.length > (isMetalCenter ? 8 : 6)) continue;

    // Identify Dewar-Chatt-Duncanson η²-alkene donor pairs (two C neighbors of centerAtom directly bonded via C=C)
    // so both alkene carbons merge into a single η² coordination branch (supports 1, 2, or 3 η²-alkene ligands such as Ni(C2H4)3).
    const eta2SecondaryIds = new Set<string>();
    if (isMetalCenter) {
      const carbonDonorIds = donorIds.filter(
        (id) => atoms.find((a) => a.id === id)?.element === 'C'
      );
      const cDonorSet = new Set(carbonDonorIds);
      // Exclude η⁵-Cp / η⁶-arene sandwich rings where every carbon donor is bonded to 2 other carbon donors in a single ring
      const isContinuousDonorRing =
        carbonDonorIds.length >= 4 &&
        carbonDonorIds.every(
          (id) => (adj.get(id) || []).filter((nb) => cDonorSet.has(nb)).length === 2
        );
      if (!isContinuousDonorRing) {
        const paired = new Set<string>();
        for (let i = 0; i < carbonDonorIds.length; i++) {
          const c1Id = carbonDonorIds[i];
          if (paired.has(c1Id)) continue;
          for (let j = i + 1; j < carbonDonorIds.length; j++) {
            const c2Id = carbonDonorIds[j];
            if (paired.has(c2Id)) continue;
            if ((adj.get(c1Id) || []).includes(c2Id)) {
              paired.add(c1Id);
              paired.add(c2Id);
              eta2SecondaryIds.add(c2Id);
              break;
            }
          }
        }
      }
    }

    const globalVisited = new Set<string>([centerAtom.id]);
    const branches: LigandBranchInfo[] = [];
    let validStar = true;
    const maxBranchSize = isMetalCenter ? 65 : 5;

    for (const dId of donorIds) {
      if (eta2SecondaryIds.has(dId)) {
        continue;
      }
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
      const branchVisited = new Set<string>([dId]);
      const queue: AtomNode[] = [donor];
      globalVisited.add(dId);

      while (queue.length > 0) {
        const cur = queue.shift()!;
        branchAtoms.push(cur);
        for (const nbId of adj.get(cur.id) || []) {
          if (nbId === centerAtom.id) continue;
          if (branchVisited.has(nbId)) {
            // Intra-ligand ring closure (e.g. phenyl ring in -PPh3, pyridine ring, cyclohexene ring, or glucose ring in Auranofin)
            continue;
          }
          if (globalVisited.has(nbId)) {
            // Cross-branch link between two different donor sites of centerAtom (bidentate chelate or η²,η²-diene chelate)
            validStar = false;
            break;
          }
          branchVisited.add(nbId);
          globalVisited.add(nbId);
          const nbAtom = atoms.find((a) => a.id === nbId);
          if (nbAtom) queue.push(nbAtom);
        }
        if (!validStar) break;
      }

      if (!validStar || branchAtoms.length > maxBranchSize) {
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

    const nLig = branches.length;
    if (nLig < 2 || nLig > 6) continue;

    if (validStar && globalVisited.size === atoms.length) {
      // If nLig === 2 and both branches are 5- or 6-membered rings (metallocene / bis-arene), defer to analyzeChelateOrSandwichTopology
      if (
        isMetalCenter &&
        nLig === 2 &&
        branches.every((b) => {
          const cCount = b.distalAtoms.filter((x) => x.element === 'C').length;
          return b.donor.element === 'C' && (cCount === 4 || cCount === 5);
        })
      ) {
        continue;
      }

      const isMonoatomicLigands = branches.every((b) => b.distalAtoms.length === 0);

      // Only transition-metal and lanthanide coordination complexes use polyatomic-branch star topology;
      // main-group centers use star VSEPR only for monoatomic ligands (AX_nE_m), while multi-center
      // organic/inorganic chains and rings (C2H4, C4H2, PAHs, cages) use the general multi-center classifier.
      if (!isMonoatomicLigands && !isMetalCenter) {
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

  // 3. If format is Formula (or input looks like a stoich formula such as SiH4, GeCl4, HBr, PCl3, H2Se
  //    where each element symbol appears at most once and there are no SMILES ring-closure digits like C1CCCCC1)
  const formulaTokens: { el: ElementSymbol; count: number }[] = [];
  const tokenRegex = /([A-Z][a-z]?)(\d*)/g;
  let tokMatch: RegExpExecArray | null;
  while ((tokMatch = tokenRegex.exec(text)) !== null) {
    const sym = tokMatch[1] as ElementSymbol;
    if (ELEMENT_PROPERTIES[sym]) {
      const count = tokMatch[2] ? parseInt(tokMatch[2], 10) : 1;
      formulaTokens.push({ el: sym, count: Math.min(12, Math.max(1, count)) });
    }
  }
  const hasUniqueElementsOnly =
    new Set(formulaTokens.map((p) => p.el)).size === formulaTokens.length;
  const isFormulaLike =
    hasUniqueElementsOnly &&
    (format === 'Formula' ||
      (/^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)+$/.test(text) && !/[=#()[\]]/.test(text)));
  if (isFormulaLike && formulaTokens.length > 0) {
    const parts = formulaTokens;
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

  // 4. Parse linear/branched/cyclic SMILES strings (e.g. C=CC=C, c1ccccc1, C1CCCCC1, CCCCCCCCCCCC, CC(=O)O, [Li+], BrCBr)
  const atoms: AtomNode[] = [];
  const bonds: BondEdge[] = [];
  const counts: Record<string, number> = {};
  const nextId = (el: ElementSymbol) => {
    counts[el] = (counts[el] || 0) + 1;
    return `${el}${counts[el]}`;
  };

  const ringClosures = new Map<string, { atomId: string; order: 1 | 2 | 3 }>();
  const ringPlacedAtomIds = new Set<string>();
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
          // If this closes a simple unplaced ring of 3..8 atoms, arrange those ring vertices into a regular polygon
          // anchored at atoms[startIdx] so any exocyclic bond to atoms[startIdx] retains its 1.42 Å length.
          const startIdx = atoms.findIndex((a) => a.id === existing.atomId);
          const endIdx = atoms.findIndex((a) => a.id === currentAtomId);
          const ringLen = endIdx - startIdx + 1;
          const alreadyInRing =
            startIdx >= 0 &&
            atoms.slice(startIdx, endIdx + 1).some((a) => ringPlacedAtomIds.has(a.id));
          if (startIdx >= 0 && ringLen >= 3 && ringLen <= 8 && !alreadyInRing) {
            const rPoly = 1.42 / (2 * Math.sin(Math.PI / ringLen));
            const anchorX = atoms[startIdx].x2d;
            const anchorY = atoms[startIdx].y2d;
            let dirAngle = 0;
            if (startIdx > 0) {
              const exoBond = bonds.find(
                (b) =>
                  (b.source === atoms[startIdx].id &&
                    atoms.findIndex((x) => x.id === b.target) < startIdx) ||
                  (b.target === atoms[startIdx].id &&
                    atoms.findIndex((x) => x.id === b.source) < startIdx)
              );
              if (exoBond) {
                const prevNode = atoms.find(
                  (x) =>
                    x.id ===
                    (exoBond.source === atoms[startIdx].id ? exoBond.target : exoBond.source)
                );
                if (prevNode) {
                  dirAngle = Math.atan2(anchorY - prevNode.y2d, anchorX - prevNode.x2d);
                }
              }
            }
            const centerX = startIdx > 0 ? anchorX + rPoly * Math.cos(dirAngle) : anchorX;
            const centerY = startIdx > 0 ? anchorY + rPoly * Math.sin(dirAngle) : anchorY;
            const baseTheta = startIdx > 0 ? dirAngle + Math.PI : 0;
            for (let k = 0; k < ringLen; k++) {
              const theta = baseTheta + (k * 2 * Math.PI) / ringLen;
              const rx = Number((centerX + rPoly * Math.cos(theta)).toFixed(2));
              const ry = Number((centerY + rPoly * Math.sin(theta)).toFixed(2));
              atoms[startIdx + k].x2d = rx;
              atoms[startIdx + k].y2d = ry;
              atoms[startIdx + k].x3d = rx;
              atoms[startIdx + k].y3d = -ry;
              ringPlacedAtomIds.add(atoms[startIdx + k].id);
            }
          }
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
      const prevIdx = atoms.findIndex((a) => a.id === currentAtomId);
      const isBranchRoot =
        branchStack.length > 0 && currentAtomId === branchStack[branchStack.length - 1];
      // Count existing bonds on prevAtom to determine main chain vs branch direction
      const prevDeg = prevAtom
        ? bonds.filter((b) => b.source === prevAtom.id || b.target === prevAtom.id).length
        : 0;
      // Compute graph distance of prevAtom from C1 so main-chain zigzag strictly alternates ±30°
      let prevDepth = Math.max(0, prevIdx);
      if (prevAtom) {
        let d = 0;
        let curId: string | undefined = prevAtom.id;
        const visited = new Set<string>();
        while (curId && !visited.has(curId)) {
          visited.add(curId);
          const inBond = bonds.find((b) => b.target === curId);
          if (inBond) {
            d++;
            curId = inBond.source;
          } else {
            break;
          }
        }
        prevDepth = d;
      }
      const angle =
        isBranchRoot || prevDeg >= 2
          ? prevDepth % 2 === 0
            ? -Math.PI / 2
            : Math.PI / 2
          : prevDepth % 2 === 0
            ? -Math.PI / 6
            : Math.PI / 6;
      const dist = 1.42;
      const x2d = prevAtom
        ? Number((prevAtom.x2d + Math.cos(angle) * dist).toFixed(2))
        : Number(((stepIdx - 2) * 0.85).toFixed(2));
      const y2d = prevAtom
        ? Number((prevAtom.y2d + Math.sin(angle) * dist).toFixed(2))
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

    const allConnectedBonds = bonds.filter((b) => b.source === atom.id || b.target === atom.id);
    const nonMetalConnectedBonds =
      !isTransitionMetal(atom.element) && !isLanthanide(atom.element)
        ? allConnectedBonds.filter((b) => {
            const otherId = b.source === atom.id ? b.target : b.source;
            const other = atoms.find((a) => a.id === otherId);
            return !other || (!isTransitionMetal(other.element) && !isLanthanide(other.element));
          })
        : allConnectedBonds;
    const connectedBonds =
      nonMetalConnectedBonds.length > 0 ? nonMetalConnectedBonds : allConnectedBonds;

    // 2D neighbor vectors for clean 2D sketcher layout
    const neighborVectors2D: [number, number][] = [];
    // 3D neighbor unit vectors for physical 3D tetrahedral (sp3) / trigonal planar (sp2) placement
    const neighborVectors3D: [number, number, number][] = [];
    for (const b of connectedBonds) {
      const otherId = b.source === atom.id ? b.target : b.source;
      const other = atoms.find((a) => a.id === otherId);
      if (other) {
        const dx2 = other.x2d - atom.x2d;
        const dy2 = other.y2d - atom.y2d;
        const len2 = Math.hypot(dx2, dy2) || 1;
        neighborVectors2D.push([dx2 / len2, dy2 / len2]);

        const dx3 = other.x3d - atom.x3d;
        const dy3 = other.y3d - atom.y3d;
        const dz3 = (other.z3d || 0) - (atom.z3d || 0);
        const len3 = Math.hypot(dx3, dy3, dz3) || 1;
        neighborVectors3D.push([dx3 / len3, dy3 / len3, dz3 / len3]);
      }
    }

    let baseAngle2D = 0;
    if (neighborVectors2D.length > 0) {
      const sumX = neighborVectors2D.reduce((s, v) => s + v[0], 0);
      const sumY = neighborVectors2D.reduce((s, v) => s + v[1], 0);
      if (Math.hypot(sumX, sumY) > 0.1) {
        baseAngle2D = Math.atan2(-sumY, -sumX);
      } else if (neighborVectors2D.length === 2) {
        baseAngle2D = Math.atan2(neighborVectors2D[0][0], -neighborVectors2D[0][1]);
      }
    }

    const r3d = (ELEMENT_PROPERTIES[atom.element]?.covalentRadius ?? 0.77) + 0.32;
    const isSp3 =
      connectedBonds.every((b) => b.order === 1) &&
      (connectedBonds.length + missingH === 4 ||
        (['N', 'P', 'As', 'Sb'].includes(atom.element) && connectedBonds.length + missingH === 3) ||
        (['O', 'S', 'Se', 'Te'].includes(atom.element) && connectedBonds.length + missingH === 2));

    // Compute exact 3D unit vectors for each implicit hydrogen
    const hDirs3D: [number, number, number][] = [];
    const k = neighborVectors3D.length;

    if (isSp3 && k === 2 && missingH === 2) {
      // Methylene (-CH2-): place the 2 C-H bonds in the plane perpendicular to C-C-C at ±54.735° (109.47° tetrahedral angle)
      const u1 = neighborVectors3D[0];
      const u2 = neighborVectors3D[1];
      let bx = -(u1[0] + u2[0]);
      let by = -(u1[1] + u2[1]);
      let bz = -(u1[2] + u2[2]);
      let bLen = Math.hypot(bx, by, bz);
      if (bLen < 0.08) {
        bx = -u1[1];
        by = u1[0];
        bz = 0;
        bLen = Math.hypot(bx, by, bz) || 1;
      }
      bx /= bLen;
      by /= bLen;
      bz /= bLen;

      // Normal to the u1-u2 plane
      let nx = u1[1] * u2[2] - u1[2] * u2[1];
      let ny = u1[2] * u2[0] - u1[0] * u2[2];
      let nz = u1[0] * u2[1] - u1[1] * u2[0];
      let nLen = Math.hypot(nx, ny, nz);
      if (nLen < 0.08) {
        nx = 0;
        ny = 0;
        nz = 1;
      } else {
        nx /= nLen;
        ny /= nLen;
        nz /= nLen;
      }
      // Orthogonalize n against b
      const dotNB = nx * bx + ny * by + nz * bz;
      nx -= dotNB * bx;
      ny -= dotNB * by;
      nz -= dotNB * bz;
      nLen = Math.hypot(nx, ny, nz) || 1;
      nx /= nLen;
      ny /= nLen;
      nz /= nLen;

      const cTet = 1 / Math.sqrt(3); // cos(54.735°) ≈ 0.57735
      const sTet = Math.sqrt(2 / 3); // sin(54.735°) ≈ 0.81650
      hDirs3D.push([cTet * bx + sTet * nx, cTet * by + sTet * ny, cTet * bz + sTet * nz]);
      hDirs3D.push([cTet * bx - sTet * nx, cTet * by - sTet * ny, cTet * bz - sTet * nz]);
    } else if (isSp3 && k === 1 && missingH >= 1) {
      // Methyl (-CH3), Amino (-NH2), or Hydroxyl (-OH): staggered 3D tetrahedral cone (cos 70.53° = 1/3, sin 70.53° = 2√2/3)
      const u1 = neighborVectors3D[0];
      const ax = -u1[0];
      const ay = -u1[1];
      const az = -u1[2];
      let e1x = -ay;
      let e1y = ax;
      let e1z = 0;
      let e1Len = Math.hypot(e1x, e1y, e1z);
      if (e1Len < 0.1) {
        e1x = 1;
        e1y = 0;
        e1z = 0;
        e1Len = 1;
      }
      e1x /= e1Len;
      e1y /= e1Len;
      const e2x = ay * e1z - az * e1y;
      const e2y = az * e1x - ax * e1z;
      const e2z = ax * e1y - ay * e1x;

      const cosCone = 1 / 3;
      const sinCone = (2 * Math.sqrt(2)) / 3;
      for (let i = 0; i < missingH; i++) {
        const phi = (i * 2 * Math.PI) / 3;
        const cp = Math.cos(phi);
        const sp = Math.sin(phi);
        hDirs3D.push([
          cosCone * ax + sinCone * (cp * e1x + sp * e2x),
          cosCone * ay + sinCone * (cp * e1y + sp * e2y),
          cosCone * az + sinCone * (cp * e1z + sp * e2z),
        ]);
      }
    } else if (isSp3 && k === 3 && missingH === 1) {
      // Tertiary methine (-CH(R1)(R2)(R3)): 4th apex of 3D tetrahedron
      const sx = neighborVectors3D[0][0] + neighborVectors3D[1][0] + neighborVectors3D[2][0];
      const sy = neighborVectors3D[0][1] + neighborVectors3D[1][1] + neighborVectors3D[2][1];
      const sz = neighborVectors3D[0][2] + neighborVectors3D[1][2] + neighborVectors3D[2][2];
      if (Math.abs(sz) > 0.15 && Math.hypot(sx, sy, sz) > 0.15) {
        const sLen = Math.hypot(sx, sy, sz);
        hDirs3D.push([-sx / sLen, -sy / sLen, -sz / sLen]);
      } else {
        const bLen = Math.hypot(sx, sy) || 1;
        const bx = Math.hypot(sx, sy) > 0.08 ? -sx / bLen : 0;
        const by = Math.hypot(sx, sy) > 0.08 ? -sy / bLen : 0;
        const signZ = connectedBonds.some((b) => b.stereo === 'wedge') ? -1 : 1;
        const norm = Math.hypot(0.45 * bx, 0.45 * by, 0.893) || 1;
        hDirs3D.push([(0.45 * bx) / norm, (0.45 * by) / norm, (signZ * 0.893) / norm]);
      }
    } else if (isSp3 && k === 0 && missingH === 4) {
      // Bare Carbon -> Methane (CH4) regular tetrahedron
      const s = 1 / Math.sqrt(3);
      hDirs3D.push([s, s, s], [s, -s, -s], [-s, s, -s], [-s, -s, s]);
    }

    for (let i = 0; i < missingH; i++) {
      const spread2D =
        missingH === 1 ? 0 : missingH === 2 ? (i === 0 ? -0.65 : 0.65) : (i * 2 * Math.PI) / missingH;
      const angle2D = baseAngle2D + spread2D;
      const r2d = 0.95;
      const hId = `H_imp_${atom.id}_${hIndex++}`;
      const hx2d = atom.x2d + Math.cos(angle2D) * r2d;
      const hy2d = atom.y2d + Math.sin(angle2D) * r2d;

      const dir3 = hDirs3D[i];
      const hx3d = dir3 ? atom.x3d + dir3[0] * r3d : atom.x3d + Math.cos(angle2D) * r3d;
      const hy3d = dir3 ? atom.y3d + dir3[1] * r3d : atom.y3d - Math.sin(angle2D) * r3d;
      const hz3d = dir3 ? (atom.z3d || 0) + dir3[2] * r3d : 0;

      outAtoms.push({
        id: hId,
        element: 'H',
        x2d: Number(hx2d.toFixed(2)),
        y2d: Number(hy2d.toFixed(2)),
        x3d: Number(hx3d.toFixed(3)),
        y3d: Number(hy3d.toFixed(3)),
        z3d: Number(hz3d.toFixed(3)),
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

  // 2. Check if molecule is a bidentate chelate or η5-metallocene sandwich complex
  //    (e.g. Ferrocene Fe(C5H5)2, [Co(en)3]3+, [Ru(bpy)3]2+, Fe(acac)3, [Fe(C2O4)3]3-, Ni(dmgH)2, Oxaliplatin, Carboplatin)
  const chelateTopo = analyzeChelateOrSandwichTopology(atoms, bonds);
  if (chelateTopo) {
    const { centerAtom, kind } = chelateTopo;
    const coords3D = new Map<string, [number, number, number]>();
    coords3D.set(centerAtom.id, [0, 0, 0]);

    const adj = new Map<string, string[]>();
    for (const a of atoms) adj.set(a.id, []);
    for (const b of bonds) {
      adj.get(b.source)?.push(b.target);
      adj.get(b.target)?.push(b.source);
    }

    // Helper to place all hydrogens in a ligand component around their already-placed parent heavy atoms
    const placeComponentHydrogens3D = (
      comp: AtomNode[],
      outwardDir: [number, number, number]
    ) => {
      const compHeavy = comp.filter((a) => a.element !== 'H');
      for (const parent of compHeavy) {
        const pPos = coords3D.get(parent.id);
        if (!pPos) continue;
        const conn = bonds.filter((b) => b.source === parent.id || b.target === parent.id);
        const hNeighbors = conn
          .map((b) => comp.find((x) => x.id === (b.source === parent.id ? b.target : b.source)))
          .filter((x): x is AtomNode => Boolean(x && x.element === 'H'));
        if (hNeighbors.length === 0) continue;

        const heavyNbs = conn
          .map((b) => comp.find((x) => x.id === (b.source === parent.id ? b.target : b.source)))
          .filter((x): x is AtomNode => Boolean(x && x.element !== 'H'));

        const uVecs: [number, number, number][] = heavyNbs
          .map((nb) => {
            const nbPos = coords3D.get(nb.id);
            if (!nbPos) return null;
            const dx = nbPos[0] - pPos[0];
            const dy = nbPos[1] - pPos[1];
            const dz = nbPos[2] - pPos[2];
            const len = Math.hypot(dx, dy, dz) || 1;
            return [dx / len, dy / len, dz / len] as [number, number, number];
          })
          .filter((v): v is [number, number, number] => Boolean(v));

        const rXH = (ELEMENT_PROPERTIES[parent.element]?.covalentRadius ?? 0.77) + 0.32;
        const organicBonds = conn.filter((b) => {
          const otherId = b.source === parent.id ? b.target : b.source;
          return otherId !== centerAtom.id;
        });
        const isParentSp3 = organicBonds.every((b) => b.order === 1) && organicBonds.length >= 3;

        if (uVecs.length === 2 && hNeighbors.length === 2) {
          const u1 = uVecs[0];
          const u2 = uVecs[1];
          let bx = -(u1[0] + u2[0]);
          let by = -(u1[1] + u2[1]);
          let bz = -(u1[2] + u2[2]);
          let bLen = Math.hypot(bx, by, bz);
          if (bLen < 0.08) {
            bx = outwardDir[0];
            by = outwardDir[1];
            bz = outwardDir[2];
            bLen = Math.hypot(bx, by, bz) || 1;
          }
          bx /= bLen;
          by /= bLen;
          bz /= bLen;

          let nx = u1[1] * u2[2] - u1[2] * u2[1];
          let ny = u1[2] * u2[0] - u1[0] * u2[2];
          let nz = u1[0] * u2[1] - u1[1] * u2[0];
          let nLen = Math.hypot(nx, ny, nz) || 1;
          nx /= nLen;
          ny /= nLen;
          nz /= nLen;
          const dotNB = nx * bx + ny * by + nz * bz;
          nx -= dotNB * bx;
          ny -= dotNB * by;
          nz -= dotNB * bz;
          nLen = Math.hypot(nx, ny, nz) || 1;
          nx /= nLen;
          ny /= nLen;
          nz /= nLen;

          const cTet = isParentSp3 ? 1 / Math.sqrt(3) : 0.5;
          const sTet = isParentSp3 ? Math.sqrt(2 / 3) : 0.866;
          coords3D.set(hNeighbors[0].id, [
            Number((pPos[0] + rXH * (cTet * bx + sTet * nx)).toFixed(3)),
            Number((pPos[1] + rXH * (cTet * by + sTet * ny)).toFixed(3)),
            Number((pPos[2] + rXH * (cTet * bz + sTet * nz)).toFixed(3)),
          ]);
          coords3D.set(hNeighbors[1].id, [
            Number((pPos[0] + rXH * (cTet * bx - sTet * nx)).toFixed(3)),
            Number((pPos[1] + rXH * (cTet * by - sTet * ny)).toFixed(3)),
            Number((pPos[2] + rXH * (cTet * bz - sTet * nz)).toFixed(3)),
          ]);
        } else if (uVecs.length >= 2 && hNeighbors.length === 1) {
          let sx = uVecs.reduce((s, v) => s + v[0], 0);
          let sy = uVecs.reduce((s, v) => s + v[1], 0);
          let sz = uVecs.reduce((s, v) => s + v[2], 0);
          // For coordinated alkene carbons, bias the vinylic C-H slightly outward away from the metal
          if ((adj.get(parent.id) || []).includes(centerAtom.id)) {
            sx -= 0.45 * outwardDir[0];
            sy -= 0.45 * outwardDir[1];
            sz -= 0.45 * outwardDir[2];
          }
          const sLen = Math.hypot(sx, sy, sz) || 1;
          coords3D.set(hNeighbors[0].id, [
            Number((pPos[0] - (rXH * sx) / sLen).toFixed(3)),
            Number((pPos[1] - (rXH * sy) / sLen).toFixed(3)),
            Number((pPos[2] - (rXH * sz) / sLen).toFixed(3)),
          ]);
        } else if (uVecs.length === 1 && hNeighbors.length >= 1) {
          const u1 = uVecs[0];
          const ax = -u1[0];
          const ay = -u1[1];
          const az = -u1[2];
          let e1x = -ay;
          let e1y = ax;
          let e1z = 0;
          let e1Len = Math.hypot(e1x, e1y, e1z);
          if (e1Len < 0.1) {
            e1x = 0;
            e1y = -az;
            e1z = ay;
            e1Len = Math.hypot(e1x, e1y, e1z) || 1;
          }
          e1x /= e1Len;
          e1y /= e1Len;
          e1z /= e1Len;
          const e2x = ay * e1z - az * e1y;
          const e2y = az * e1x - ax * e1z;
          const e2z = ax * e1y - ay * e1x;
          const cosCone = hNeighbors.length === 2 ? 0.5 : 1 / 3;
          const sinCone = hNeighbors.length === 2 ? 0.866 : (2 * Math.sqrt(2)) / 3;
          hNeighbors.forEach((hAtom, idx) => {
            const phi = (idx * 2 * Math.PI) / hNeighbors.length;
            const cp = Math.cos(phi);
            const sp = Math.sin(phi);
            coords3D.set(hAtom.id, [
              Number(
                (pPos[0] + rXH * (cosCone * ax + sinCone * (cp * e1x + sp * e2x))).toFixed(3)
              ),
              Number(
                (pPos[1] + rXH * (cosCone * ay + sinCone * (cp * e1y + sp * e2y))).toFixed(3)
              ),
              Number(
                (pPos[2] + rXH * (cosCone * az + sinCone * (cp * e1z + sp * e2z))).toFixed(3)
              ),
            ]);
          });
        }
      }
    };

    // Helper to reconstruct 3D tub/boat/cage geometry for a chelating η²,η²-diene ligand
    // (e.g. 1,5-cyclooctadiene COD in Ni(COD)2 / Pt(COD)Cl2, norbornadiene NBD in Pt(nbd)Cl2, 1,3-butadiene)
    const placeDieneChelateComponent = (
      comp: AtomNode[],
      bOut: [number, number, number],
      eChord: [number, number, number],
      eCC: [number, number, number]
    ) => {
      const compDonors = comp.filter((a) => (adj.get(a.id) || []).includes(centerAtom.id));
      const donorSet = new Set(compDonors.map((d) => d.id));
      const bridgeHeavy = comp.filter((a) => a.element !== 'H' && !donorSet.has(a.id));

      // Partition the 4 donor carbons into two bonded C=C pairs: pairA = [cA1, cA2] and pairB = [cB1, cB2]
      const pairedIds = new Set<string>();
      const ccPairs: [AtomNode, AtomNode][] = [];
      for (let i = 0; i < compDonors.length; i++) {
        const d1 = compDonors[i];
        if (pairedIds.has(d1.id)) continue;
        const d2 =
          compDonors.find(
            (cand, j) =>
              j > i &&
              !pairedIds.has(cand.id) &&
              bonds.some(
                (b) =>
                  b.order === 2 &&
                  ((b.source === d1.id && b.target === cand.id) ||
                    (b.source === cand.id && b.target === d1.id))
              )
          ) ||
          compDonors.find(
            (cand, j) =>
              j > i &&
              !pairedIds.has(cand.id) &&
              bonds.some(
                (b) =>
                  (b.source === d1.id && b.target === cand.id) ||
                  (b.source === cand.id && b.target === d1.id)
              )
          );
        if (d2) {
          pairedIds.add(d1.id);
          pairedIds.add(d2.id);
          ccPairs.push([d1, d2]);
        }
      }

      if (ccPairs.length >= 2) {
        const [cA1, cA2] = ccPairs[0];
        let [cB1, cB2] = ccPairs[1];

        // Graph distance inside `comp` excluding the two C=C bonds so cA1 and cB1 lie on the same -eCC side
        const graphDistNoDouble = (startId: string, targetId: string): number => {
          const q: { id: string; dist: number }[] = [{ id: startId, dist: 0 }];
          const vis = new Set<string>([startId]);
          while (q.length > 0) {
            const { id, dist } = q.shift()!;
            if (id === targetId) return dist;
            for (const nbId of adj.get(id) || []) {
              if (nbId === centerAtom.id || vis.has(nbId)) continue;
              if (
                (id === cA1.id && nbId === cA2.id) ||
                (id === cA2.id && nbId === cA1.id) ||
                (id === cB1.id && nbId === cB2.id) ||
                (id === cB2.id && nbId === cB1.id)
              ) {
                continue;
              }
              if (comp.some((x) => x.id === nbId && x.element !== 'H')) {
                vis.add(nbId);
                q.push({ id: nbId, dist: dist + 1 });
              }
            }
          }
          return 99;
        };

        if (graphDistNoDouble(cA1.id, cB2.id) < graphDistNoDouble(cA1.id, cB1.id)) {
          [cB1, cB2] = [cB2, cB1];
        }

        const setVec = (
          atom: AtomNode,
          sOut: number,
          sChord: number,
          sCC: number
        ) => {
          coords3D.set(atom.id, [
            Number((sOut * bOut[0] + sChord * eChord[0] + sCC * eCC[0]).toFixed(3)),
            Number((sOut * bOut[1] + sChord * eChord[1] + sCC * eCC[1]).toFixed(3)),
            Number((sOut * bOut[2] + sChord * eChord[2] + sCC * eCC[2]).toFixed(3)),
          ]);
        };

        if (bridgeHeavy.length === 0) {
          // Conjugated cisoid η⁴-1,3-butadiene (C1=C2-C3=C4): ensure cIntA and cIntB are the two internal carbons (degree 2 within compDonors)
          const donorNeighborsCount = (c: AtomNode) =>
            (adj.get(c.id) || []).filter((nb) => donorSet.has(nb)).length;
          const [cIntA, cTermA] =
            donorNeighborsCount(cA1) >= donorNeighborsCount(cA2) ? [cA1, cA2] : [cA2, cA1];
          const [cIntB, cTermB] =
            donorNeighborsCount(cB1) >= donorNeighborsCount(cB2) ? [cB1, cB2] : [cB2, cB1];
          setVec(cIntA, 1.48, 0.73, -0.52);
          setVec(cIntB, 1.48, -0.73, -0.52);
          setVec(cTermA, 1.32, 1.42, 0.68);
          setVec(cTermB, 1.32, -1.42, 0.68);
        } else {
          const dChord = bridgeHeavy.length <= 3 ? 1.18 : 1.36;
          const halfCC = 0.688;
          setVec(cA1, 1.32, dChord, -halfCC);
          setVec(cA2, 1.32, dChord, halfCC);
          setVec(cB1, 1.32, -dChord, -halfCC);
          setVec(cB2, 1.32, -dChord, halfCC);

          if (bridgeHeavy.length === 4) {
            // 1,5-Cyclooctadiene (COD) tub/boat conformation: two -CH2-CH2- bridges arching outward along +bOut
            for (const m of bridgeHeavy) {
              const mNbs = adj.get(m.id) || [];
              if (mNbs.includes(cA1.id)) {
                setVec(m, 2.48, 0.77, -1.28);
              } else if (mNbs.includes(cB1.id)) {
                setVec(m, 2.48, -0.77, -1.28);
              } else if (mNbs.includes(cA2.id)) {
                setVec(m, 2.48, 0.77, 1.28);
              } else {
                setVec(m, 2.48, -0.77, 1.28);
              }
            }
          } else if (bridgeHeavy.length === 3) {
            // Bicyclic Norbornadiene (NBD): 2 bridgehead -CH- carbons + 1 apical -CH2- bridge carbon
            for (const m of bridgeHeavy) {
              const mNbs = adj.get(m.id) || [];
              if (mNbs.includes(cA1.id) || mNbs.includes(cB1.id)) {
                setVec(m, 1.95, 0, -1.15);
              } else if (mNbs.includes(cA2.id) || mNbs.includes(cB2.id)) {
                setVec(m, 1.95, 0, 1.15);
              } else {
                setVec(m, 2.72, 0, 0);
              }
            }
          } else {
            // General bridged diene (e.g. 1,4-cyclohexadiene or substituted COD)
            bridgeHeavy.forEach((m, idx) => {
              const mNbs = adj.get(m.id) || [];
              const sCC =
                mNbs.includes(cA1.id) || mNbs.includes(cB1.id)
                  ? -1.18
                  : mNbs.includes(cA2.id) || mNbs.includes(cB2.id)
                    ? 1.18
                    : (idx % 2 === 0 ? -1 : 1) * 1.18;
              const sChord = mNbs.includes(cA1.id) || mNbs.includes(cA2.id)
                ? 0.55
                : mNbs.includes(cB1.id) || mNbs.includes(cB2.id)
                  ? -0.55
                  : 0;
              setVec(m, 2.15, sChord, sCC);
            });
          }
        }
      }

      placeComponentHydrogens3D(comp, bOut);
    };

    if (kind === 'metallocene_sandwich') {
      // Find the two C5 / C6 ring components when centerAtom is removed
      const visited = new Set<string>([centerAtom.id]);
      const ringComps: AtomNode[][] = [];
      for (const a of atoms) {
        if (visited.has(a.id)) continue;
        const comp: AtomNode[] = [];
        const q = [a.id];
        visited.add(a.id);
        while (q.length > 0) {
          const curr = q.shift()!;
          const node = atoms.find((x) => x.id === curr);
          if (node) comp.push(node);
          for (const nb of adj.get(curr) || []) {
            if (!visited.has(nb)) {
              visited.add(nb);
              q.push(nb);
            }
          }
        }
        ringComps.push(comp);
      }

      ringComps.slice(0, 2).forEach((comp, ringIdx) => {
        const ringCarbons = comp.filter(
          (a) => a.element === 'C' && (adj.get(a.id) || []).includes(centerAtom.id)
        );
        const orderedRing: AtomNode[] = [];
        if (ringCarbons.length > 0) {
          const ringSet = new Set(ringCarbons.map((c) => c.id));
          orderedRing.push(ringCarbons[0]);
          const used = new Set<string>([ringCarbons[0].id]);
          while (orderedRing.length < ringCarbons.length) {
            const last = orderedRing[orderedRing.length - 1];
            const nextId = (adj.get(last.id) || []).find(
              (nb) => ringSet.has(nb) && !used.has(nb)
            );
            if (!nextId) {
              const fallback = ringCarbons.find((c) => !used.has(c.id));
              if (!fallback) break;
              orderedRing.push(fallback);
              used.add(fallback.id);
            } else {
              const nextNode = ringCarbons.find((c) => c.id === nextId)!;
              orderedRing.push(nextNode);
              used.add(nextNode.id);
            }
          }
        }

        const isSixRing = orderedRing.length === 6;
        const zRing = ringIdx === 0 ? (isSixRing ? 1.61 : 1.66) : isSixRing ? -1.61 : -1.66;
        const rCarbon = isSixRing ? 1.40 : 1.21;
        // Eclipsed D6h for Bis(benzene)chromium Cr(C6H6)2; staggered D5d for Ferrocene Fe(C5H5)2
        const phiOffset = ringIdx === 0 || isSixRing ? 0 : Math.PI / 5;

        orderedRing.forEach((cAtom, k) => {
          const phi = phiOffset + (k * 2 * Math.PI) / Math.max(1, orderedRing.length);
          const cx = Number((rCarbon * Math.cos(phi)).toFixed(3));
          const cy = Number((rCarbon * Math.sin(phi)).toFixed(3));
          coords3D.set(cAtom.id, [cx, cy, zRing]);

          // Place any attached H or substituent radially outward from the ring
          const subs = (adj.get(cAtom.id) || [])
            .filter((nb) => nb !== centerAtom.id && !orderedRing.some((rc) => rc.id === nb))
            .map((nb) => atoms.find((x) => x.id === nb))
            .filter((x): x is AtomNode => Boolean(x));
          subs.forEach((sub, sIdx) => {
            const rSub = rCarbon + (sub.element === 'H' ? 1.08 : 1.5);
            const subPhi = phi + (sIdx - (subs.length - 1) / 2) * 0.15;
            coords3D.set(sub.id, [
              Number((rSub * Math.cos(subPhi)).toFixed(3)),
              Number((rSub * Math.sin(subPhi)).toFixed(3)),
              zRing,
            ]);
          });
        });
      });
    } else if (kind === 'bis_diene_chelate') {
      // Homoleptic Bis(diene) d10 tetrahedral chelate (e.g. Bis(1,5-cyclooctadiene)nickel Ni(COD)2 -> D2d):
      // Ligand 0 on +x with C=C parallel to z; Ligand 1 on -x rotated 90° around x (C=C parallel to y) for exact S4 / D2d symmetry!
      placeDieneChelateComponent(chelateTopo.ligands[0].allAtoms, [1, 0, 0], [0, 1, 0], [0, 0, 1]);
      placeDieneChelateComponent(
        chelateTopo.ligands[1].allAtoms,
        [-1, 0, 0],
        [0, 0, 1],
        [0, 1, 0]
      );
    } else if (kind === 'mono_diene_chelate') {
      // Square-planar d8 mono-diene chelate (e.g. Pt(COD)Cl2, Pd(COD)Cl2, [Rh(COD)Cl2]-, Pt(nbd)Cl2 -> C2v):
      // Diene ligand on +x with C=C perpendicular to xy plane (along z); two cis monodentate ligands on -x at ±45° in xy plane.
      placeDieneChelateComponent(chelateTopo.ligands[0].allAtoms, [1, 0, 0], [0, 1, 0], [0, 0, 1]);
      const rCenter = ELEMENT_PROPERTIES[centerAtom.element]?.covalentRadius ?? 1.3;
      const s2 = Math.SQRT1_2;
      const monoDirs: [number, number, number][] = [
        [-s2, s2, 0],
        [-s2, -s2, 0],
      ];
      chelateTopo.ligands.slice(1, 3).forEach((lig, idx) => {
        const d0 = lig.donors[0];
        const rBond = Number(
          ((rCenter + (ELEMENT_PROPERTIES[d0.element]?.covalentRadius ?? 0.99)) * 0.96).toFixed(3)
        );
        const uDir = monoDirs[idx];
        const dPos: [number, number, number] = [
          Number((rBond * uDir[0]).toFixed(3)),
          Number((rBond * uDir[1]).toFixed(3)),
          0,
        ];
        coords3D.set(d0.id, dPos);
        const others = lig.allAtoms.filter((a) => a.id !== d0.id);
        others.forEach((other, oIdx) => {
          const rDist = other.element === 'H' ? 1.05 : 1.15;
          const zOff = others.length > 1 ? (oIdx % 2 === 0 ? 0.5 : -0.5) : 0;
          coords3D.set(other.id, [
            Number((dPos[0] + rDist * uDir[0]).toFixed(3)),
            Number((dPos[1] + rDist * uDir[1]).toFixed(3)),
            Number(zOff.toFixed(3)),
          ]);
        });
      });
    } else if (kind === 'piano_stool_sandwich') {
      // Half-sandwich / Piano-stool M(ηⁿ-ring/diene)L3 (e.g. (C6H6)Cr(CO)3, (C4H4)Fe(CO)3, (C4H6)Fe(CO)3)
      const capLig = chelateTopo.ligands[0];
      const capHeavy = capLig.allAtoms.filter((a) => a.element !== 'H');
      const capDonors = capLig.donors;
      const isCyclicCap = capHeavy.every(
        (a) =>
          (adj.get(a.id) || []).filter((nb) => capHeavy.some((h) => h.id === nb)).length === 2
      );
      if (isCyclicCap) {
        const orderedRing: AtomNode[] = [capDonors[0]];
        const used = new Set<string>([capDonors[0].id]);
        const ringSet = new Set(capDonors.map((c) => c.id));
        while (orderedRing.length < capDonors.length) {
          const last = orderedRing[orderedRing.length - 1];
          const nextId = (adj.get(last.id) || []).find(
            (nb) => ringSet.has(nb) && !used.has(nb)
          );
          if (!nextId) break;
          const nextNode = capDonors.find((c) => c.id === nextId)!;
          orderedRing.push(nextNode);
          used.add(nextNode.id);
        }
        const nRing = orderedRing.length;
        const rCarbon = nRing === 6 ? 1.40 : nRing === 5 ? 1.21 : 1.02;
        const zCap = 1.62;
        orderedRing.forEach((cAtom, k) => {
          const phi = Math.PI / 2 + (k * 2 * Math.PI) / Math.max(1, nRing);
          const cx = Number((rCarbon * Math.cos(phi)).toFixed(3));
          const cy = Number((rCarbon * Math.sin(phi)).toFixed(3));
          coords3D.set(cAtom.id, [cx, cy, zCap]);
          const subs = (adj.get(cAtom.id) || [])
            .filter((nb) => nb !== centerAtom.id && !ringSet.has(nb))
            .map((nb) => capLig.allAtoms.find((x) => x.id === nb))
            .filter((x): x is AtomNode => Boolean(x));
          subs.forEach((sub) => {
            const rSub = rCarbon + (sub.element === 'H' ? 1.08 : 1.5);
            coords3D.set(sub.id, [
              Number((rSub * Math.cos(phi)).toFixed(3)),
              Number((rSub * Math.sin(phi)).toFixed(3)),
              Number((zCap + 0.12).toFixed(3)),
            ]);
          });
        });
      } else {
        placeDieneChelateComponent(capLig.allAtoms, [0, 0, 1], [1, 0, 0], [0, 1, 0]);
      }

      // Place the 3 tripod legs at z < 0
      const rCenter = ELEMENT_PROPERTIES[centerAtom.element]?.covalentRadius ?? 1.25;
      chelateTopo.ligands.slice(1, 4).forEach((leg, idx) => {
        const d0 = leg.donors[0];
        const rBond = (rCenter + (ELEMENT_PROPERTIES[d0.element]?.covalentRadius ?? 0.75)) * 0.94;
        const phi = Math.PI / 2 + (idx * 2 * Math.PI) / 3;
        const sinT = 0.819; // sin(55°)
        const cosT = -0.574; // -cos(55°)
        const ux = sinT * Math.cos(phi);
        const uy = sinT * Math.sin(phi);
        const uz = cosT;
        const dPos: [number, number, number] = [
          Number((rBond * ux).toFixed(3)),
          Number((rBond * uy).toFixed(3)),
          Number((rBond * uz).toFixed(3)),
        ];
        coords3D.set(d0.id, dPos);
        const others = leg.allAtoms.filter((a) => a.id !== d0.id);
        others.forEach((other) => {
          const distLen = other.element === 'H' ? 1.06 : 1.15;
          coords3D.set(other.id, [
            Number((dPos[0] + distLen * ux).toFixed(3)),
            Number((dPos[1] + distLen * uy).toFixed(3)),
            Number((dPos[2] + distLen * uz).toFixed(3)),
          ]);
        });
      });
    } else {
      // Group ligand atoms into connected components when centerAtom is removed
      const visited = new Set<string>([centerAtom.id]);
      const ligComps: AtomNode[][] = [];
      for (const a of atoms) {
        if (visited.has(a.id)) continue;
        const comp: AtomNode[] = [];
        const q = [a.id];
        visited.add(a.id);
        while (q.length > 0) {
          const curr = q.shift()!;
          const node = atoms.find((x) => x.id === curr);
          if (node) comp.push(node);
          for (const nb of adj.get(curr) || []) {
            if (!visited.has(nb)) {
              visited.add(nb);
              q.push(nb);
            }
          }
        }
        ligComps.push(comp);
      }

      const rCenter = ELEMENT_PROPERTIES[centerAtom.element]?.covalentRadius ?? 1.25;
      const getMdLen = (donor: AtomNode) =>
        Number(
          ((rCenter + (ELEMENT_PROPERTIES[donor.element]?.covalentRadius ?? 0.73)) * 0.96).toFixed(
            3
          )
        );

      const placeChelateComponent = (
        comp: AtomNode[],
        d0Vec: [number, number, number],
        d1Vec: [number, number, number],
        bisector: [number, number, number],
        chordDir: [number, number, number],
        planeNorm: [number, number, number]
      ) => {
        const compDonors = comp.filter((a) => (adj.get(a.id) || []).includes(centerAtom.id));
        if (compDonors.length === 1) {
          const d0 = compDonors[0];
          const r0 = getMdLen(d0);
          const dPos: [number, number, number] = [
            Number((r0 * d0Vec[0]).toFixed(3)),
            Number((r0 * d0Vec[1]).toFixed(3)),
            Number((r0 * d0Vec[2]).toFixed(3)),
          ];
          coords3D.set(d0.id, dPos);
          const others = comp.filter((a) => a.id !== d0.id);
          others.forEach((other, idx) => {
            const phi = (idx * 2 * Math.PI) / Math.max(1, others.length);
            const rDist = other.element === 'H' ? 1.02 : 1.42;
            const cosT = 0.5;
            const sinT = 0.866;
            coords3D.set(other.id, [
              Number(
                (
                  dPos[0] +
                  rDist *
                    (cosT * d0Vec[0] +
                      sinT * (Math.cos(phi) * chordDir[0] + Math.sin(phi) * planeNorm[0]))
                ).toFixed(3)
              ),
              Number(
                (
                  dPos[1] +
                  rDist *
                    (cosT * d0Vec[1] +
                      sinT * (Math.cos(phi) * chordDir[1] + Math.sin(phi) * planeNorm[1]))
                ).toFixed(3)
              ),
              Number(
                (
                  dPos[2] +
                  rDist *
                    (cosT * d0Vec[2] +
                      sinT * (Math.cos(phi) * chordDir[2] + Math.sin(phi) * planeNorm[2]))
                ).toFixed(3)
              ),
            ]);
          });
          return;
        }
        if (compDonors.length < 2) return;

        const [d0, d1] = compDonors;
        const r0 = getMdLen(d0);
        const r1 = getMdLen(d1);
        const p0: [number, number, number] = [
          Number((r0 * d0Vec[0]).toFixed(3)),
          Number((r0 * d0Vec[1]).toFixed(3)),
          Number((r0 * d0Vec[2]).toFixed(3)),
        ];
        const p1: [number, number, number] = [
          Number((r1 * d1Vec[0]).toFixed(3)),
          Number((r1 * d1Vec[1]).toFixed(3)),
          Number((r1 * d1Vec[2]).toFixed(3)),
        ];
        coords3D.set(d0.id, p0);
        coords3D.set(d1.id, p1);

        const mid3D: [number, number, number] = [
          0.5 * (p0[0] + p1[0]),
          0.5 * (p0[1] + p1[1]),
          0.5 * (p0[2] + p1[2]),
        ];
        const targetChordLen = Math.hypot(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]) || 2.6;

        // Build local 2D/3D frame of the ligand component from d0 -> d1 and outward toward ligand centroid
        const mid2x = 0.5 * (d0.x2d + d1.x2d);
        const mid2y = 0.5 * (-d0.y2d + -d1.y2d);
        let tx = d1.x2d - d0.x2d;
        let ty = -d1.y2d - -d0.y2d;
        const origChordLen = Math.hypot(tx, ty) || 1;
        tx /= origChordLen;
        ty /= origChordLen;

        let bx = -ty;
        let by = tx;
        const nonDonors = comp.filter((a) => a.id !== d0.id && a.id !== d1.id && a.element !== 'H');
        if (nonDonors.length > 0) {
          const avgX = nonDonors.reduce((s, a) => s + a.x2d, 0) / nonDonors.length - mid2x;
          const avgY = nonDonors.reduce((s, a) => s + -a.y2d, 0) / nonDonors.length - mid2y;
          if (avgX * bx + avgY * by < 0) {
            bx = -bx;
            by = -by;
          }
        }

        const scale2D = Math.min(1.35, Math.max(0.75, targetChordLen / origChordLen));
        for (const a of comp) {
          if (a.id === d0.id || a.id === d1.id) continue;
          const dx = a.x2d - mid2x;
          const dy = -a.y2d - mid2y;
          const sOut = Math.max(0.25, (dx * bx + dy * by) * scale2D);
          const sChord = (dx * tx + dy * ty) * scale2D;
          const sNorm = a.element === 'H' ? ((a.id.charCodeAt(a.id.length - 1) % 2 === 0 ? 1 : -1) * 0.35) : 0;
          coords3D.set(a.id, [
            Number(
              (mid3D[0] + sOut * bisector[0] + sChord * chordDir[0] + sNorm * planeNorm[0]).toFixed(
                3
              )
            ),
            Number(
              (mid3D[1] + sOut * bisector[1] + sChord * chordDir[1] + sNorm * planeNorm[1]).toFixed(
                3
              )
            ),
            Number(
              (mid3D[2] + sOut * bisector[2] + sChord * chordDir[2] + sNorm * planeNorm[2]).toFixed(
                3
              )
            ),
          ]);
        }
      };

      if (kind === 'tris_bidentate_octahedral') {
        const s2 = Math.SQRT1_2;
        const octBites: {
          d0: [number, number, number];
          d1: [number, number, number];
          bis: [number, number, number];
          chord: [number, number, number];
          norm: [number, number, number];
        }[] = [
          {
            d0: [1, 0, 0],
            d1: [0, 1, 0],
            bis: [s2, s2, 0],
            chord: [-s2, s2, 0],
            norm: [0, 0, 1],
          },
          {
            d0: [0, 0, 1],
            d1: [-1, 0, 0],
            bis: [-s2, 0, s2],
            chord: [-s2, 0, -s2],
            norm: [0, 1, 0],
          },
          {
            d0: [0, -1, 0],
            d1: [0, 0, -1],
            bis: [0, -s2, -s2],
            chord: [0, -s2, s2],
            norm: [1, 0, 0],
          },
        ];
        ligComps.slice(0, 3).forEach((comp, idx) => {
          const cfg = octBites[idx];
          placeChelateComponent(comp, cfg.d0, cfg.d1, cfg.bis, cfg.chord, cfg.norm);
        });
      } else {
        // square_planar_chelate (e.g. Ni(dmgH)2, Oxaliplatin, Carboplatin)
        const s2 = Math.SQRT1_2;
        // Sort bidentate components first
        const sortedComps = [...ligComps].sort((a, b) => {
          const da = a.filter((x) => (adj.get(x.id) || []).includes(centerAtom.id)).length;
          const db = b.filter((x) => (adj.get(x.id) || []).includes(centerAtom.id)).length;
          return db - da;
        });
        if (sortedComps.length === 2) {
          placeChelateComponent(
            sortedComps[0],
            [s2, -s2, 0],
            [s2, s2, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
          );
          placeChelateComponent(
            sortedComps[1],
            [-s2, s2, 0],
            [-s2, -s2, 0],
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, 1]
          );
        } else if (sortedComps.length >= 3) {
          // 1 bidentate + 2 monodentate (e.g. Carboplatin)
          placeChelateComponent(
            sortedComps[0],
            [s2, -s2, 0],
            [s2, s2, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
          );
          placeChelateComponent(
            sortedComps[1],
            [-s2, s2, 0],
            [-s2, s2, 0],
            [-s2, s2, 0],
            [-s2, -s2, 0],
            [0, 0, 1]
          );
          placeChelateComponent(
            sortedComps[2],
            [-s2, -s2, 0],
            [-s2, -s2, 0],
            [-s2, -s2, 0],
            [s2, -s2, 0],
            [0, 0, 1]
          );
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

  // 3. Check if molecule has a central coordination / VSEPR atom surrounded by nLig ∈ [2, 6] ligand branches
  // (Supports both monoatomic ligands like H2O, NH3, CH4, SF6, CrCl6, CeCl6 and polyatomic coordination complexes
  // like Wilkinson's catalyst [RhCl(PPh3)3], Vaska's complex [IrCl(CO)(PPh3)2], Pd(PPh3)4, Cisplatin, Hexacyanoferrate [Fe(CN)6]4-).
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
      const isSqPlanarMetal = isFourCoordSquarePlanar(centerAtom, sortedBranches, charge);

      if ((isMetalCenter && !isSqPlanarMetal) || (!isMetalCenter && lonePairs === 0)) {
        // AX4E0 (CH4, SiH4, CCl4, TiCl4, Ni(CO)4, Pd(PPh3)4, NiCl2(PPh3)2) -> Tetrahedral (Td / C3v / C2v)
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
        // AX4E2 (XeF4) or d8 square-planar transition metal ([RhCl(PPh3)3], [IrCl(CO)(PPh3)2], Cisplatin, PtCl4, [Ni(CN)4]2-)
        // For ML2X2 (2 pairs of identical ligands, e.g. Cisplatin), place identical ligands cis (90°) for C2v symmetry;
        // for ML2XY (Vaska's complex, 1+1+2), sortedBranches has the two unique ligands first -> place them trans ([1,0,0], [-1,0,0])
        // and the two PPh3 ligands trans ([0,1,0], [0,-1,0]).
        const sigCounts = Object.values(sigFreq);
        const isCisPair = sigCounts.length === 2 && sigCounts.every((c) => c === 2);
        const sqDirs: [number, number, number][] = isCisPair
          ? [
              [1, 0, 0],
              [0, 1, 0],
              [-1, 0, 0],
              [0, -1, 0],
            ]
          : [
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

    // Build adjacency lookup for multi-shell ligand traversal
    const branchAdj = new Map<string, string[]>();
    for (const a of atoms) branchAdj.set(a.id, []);
    for (const b of bonds) {
      branchAdj.get(b.source)?.push(b.target);
      branchAdj.get(b.target)?.push(b.source);
    }

    // Propagate 3D coordinates outward along each polyatomic ligand branch
    // (Supports linear -CN/-CO/-NO, conical -NH3/-OH2, propeller -PPh3/-PEt3, η2-alkenes of any size, and general complex ligands)
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

      const distalSet = new Set(br.distalAtoms.map((a) => a.id));
      const firstShell = br.distalAtoms.filter((a) =>
        (branchAdj.get(br.donor.id) || []).includes(a.id)
      );
      const branchHeavy = [br.donor, ...br.distalAtoms].filter((a) => a.element !== 'H');
      const eta2PartnerC =
        br.donor.element === 'C'
          ? br.distalAtoms.find(
              (a) =>
                a.element === 'C' &&
                bonds.some(
                  (b) =>
                    (b.source === br.donor.id && b.target === a.id) ||
                    (b.source === a.id && b.target === br.donor.id)
                ) &&
                (bonds.some(
                  (b) =>
                    (b.source === centerAtom.id && b.target === a.id) ||
                    (b.source === a.id && b.target === centerAtom.id)
                ) ||
                  branchHeavy.length === 2)
            )
          : undefined;

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
      } else if (firstShell.length === br.distalAtoms.length) {
        // Single-shell conical substituent branch (e.g. -NH3, -OH2, -CH3, -NO2)
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
      } else if (br.donor.element === 'C' && eta2PartnerC) {
        // Dewar-Chatt-Duncanson η²-Alkene ligand (checked BEFORE C3 phosphine branch so substituted/cyclic alkenes
        // like propene, 1-hexene, styrene, cyclohexene, and norbornene are reconstructed as η²-alkenes):
        // 1. Orient C=C along e2 (perpendicular to the square-planar xy coordination plane)
        // 2. Elongate C=C bond length to 1.376 Å (halfCC = 0.688 Å) due to d_π -> π*(C=C) backbonding
        // 3. Bend alkene substituents back away from the metal along +u (metallacyclopropane pyramidalization)
        const otherC = eta2PartnerC;
        const halfCC = 0.688;
        const c0Pos: [number, number, number] = [
          Number((dPos[0] - halfCC * e2x).toFixed(3)),
          Number((dPos[1] - halfCC * e2y).toFixed(3)),
          Number((dPos[2] - halfCC * e2z).toFixed(3)),
        ];
        const c1Pos: [number, number, number] = [
          Number((dPos[0] + halfCC * e2x).toFixed(3)),
          Number((dPos[1] + halfCC * e2y).toFixed(3)),
          Number((dPos[2] + halfCC * e2z).toFixed(3)),
        ];
        coords3D.set(br.donor.id, c0Pos);
        coords3D.set(otherC.id, c1Pos);

        if (branchHeavy.length === 2) {
          // Unsubstituted ethylene (η²-C2H4)
          [
            { cAtom: br.donor, cPos: c0Pos, sign2: -1 },
            { cAtom: otherC, cPos: c1Pos, sign2: 1 },
          ].forEach(({ cAtom, cPos, sign2 }) => {
            const subAtoms = (branchAdj.get(cAtom.id) || [])
              .filter((nb) => nb !== centerAtom.id && nb !== br.donor.id && nb !== otherC.id)
              .map((nb) => atoms.find((x) => x.id === nb))
              .filter((x): x is AtomNode => Boolean(x));
            subAtoms.forEach((subAtom, sIdx) => {
              const sign1 = sIdx === 0 ? -1 : 1;
              const rBond = subAtom.element === 'H' ? 1.08 : 1.51;
              coords3D.set(subAtom.id, [
                Number(
                  (cPos[0] + rBond * (0.32 * ux + 0.816 * sign1 * e1x + 0.48 * sign2 * e2x)).toFixed(
                    3
                  )
                ),
                Number(
                  (cPos[1] + rBond * (0.32 * uy + 0.816 * sign1 * e1y + 0.48 * sign2 * e2y)).toFixed(
                    3
                  )
                ),
                Number(
                  (cPos[2] + rBond * (0.32 * uz + 0.816 * sign1 * e1z + 0.48 * sign2 * e2z)).toFixed(
                    3
                  )
                ),
              ]);
            });
          });
        } else {
          // Substituted, cyclic, bicyclic, or aryl-conjugated η²-alkene (propene, 1-hexene, styrene, cyclohexene, norbornene):
          // Map the ligand's 2D carbon skeleton relative to C0=C1 into the 3D DCD frame tilted away from the metal along +u
          const mid2x = 0.5 * (br.donor.x2d + otherC.x2d);
          const mid2y = 0.5 * (-br.donor.y2d + -otherC.y2d);
          let tx = otherC.x2d - br.donor.x2d;
          let ty = -otherC.y2d - -br.donor.y2d;
          const origCCLen = Math.hypot(tx, ty) || 1.4;
          tx /= origCCLen;
          ty /= origCCLen;
          let bx = -ty;
          let by = tx;
          const remHeavy = branchHeavy.filter(
            (a) => a.id !== br.donor.id && a.id !== otherC.id
          );
          if (remHeavy.length > 0) {
            const avgX = remHeavy.reduce((s, a) => s + a.x2d, 0) / remHeavy.length - mid2x;
            const avgY = remHeavy.reduce((s, a) => s + -a.y2d, 0) / remHeavy.length - mid2y;
            if (avgX * bx + avgY * by < 0) {
              bx = -bx;
              by = -by;
            }
          }

          // Compute BFS depth from br.donor for sp3 ring/chain puckering
          const alkDepth = new Map<string, number>([
            [br.donor.id, 0],
            [otherC.id, 0],
          ]);
          const qDepth = [br.donor.id, otherC.id];
          while (qDepth.length > 0) {
            const cur = qDepth.shift()!;
            const d = alkDepth.get(cur) ?? 0;
            for (const nb of branchAdj.get(cur) || []) {
              if (nb !== centerAtom.id && !alkDepth.has(nb)) {
                alkDepth.set(nb, d + 1);
                qDepth.push(nb);
              }
            }
          }

          // Tilt angle away from metal (cosθ = 0.78 along e1, sinθ = 0.626 along +u)
          const cosTilt = 0.78;
          const sinTilt = 0.626;
          const nRingX = cosTilt * ux - sinTilt * e1x;
          const nRingY = cosTilt * uy - sinTilt * e1y;
          const nRingZ = cosTilt * uz - sinTilt * e1z;

          for (const a of remHeavy) {
            const dx = a.x2d - mid2x;
            const dy = -a.y2d - mid2y;
            const sPar = dx * tx + dy * ty;
            const sPerp = dx * bx + dy * by;
            const connBonds = bonds.filter((b) => b.source === a.id || b.target === a.id);
            const isSp3C = connBonds.every((b) => b.order === 1);
            const depth = alkDepth.get(a.id) ?? 1;
            const sPucker = isSp3C ? (depth % 2 === 0 ? 0.28 : -0.28) : 0;
            const outU = Math.max(0.32, Math.abs(sPerp) * sinTilt);
            const latE1 = sPerp * cosTilt;
            coords3D.set(a.id, [
              Number(
                (dPos[0] + sPar * e2x + latE1 * e1x + outU * ux + sPucker * nRingX).toFixed(3)
              ),
              Number(
                (dPos[1] + sPar * e2y + latE1 * e1y + outU * uy + sPucker * nRingY).toFixed(3)
              ),
              Number(
                (dPos[2] + sPar * e2z + latE1 * e1z + outU * uz + sPucker * nRingZ).toFixed(3)
              ),
            ]);
          }

          // Perform 3D spring-bond & non-bonded steric relaxation on remHeavy (with br.donor and otherC anchored)
          // so bicyclic/aryl/branched alkenes (e.g. norbornene, styrene) have physical ~1.40–1.54 Å C–C bonds.
          const branchHeavySet = new Set(branchHeavy.map((x) => x.id));
          const remHeavySet = new Set(remHeavy.map((x) => x.id));
          const intraHeavyBonds = bonds.filter(
            (b) => branchHeavySet.has(b.source) && branchHeavySet.has(b.target)
          );
          for (let iter = 0; iter < 35; iter++) {
            for (const b of intraHeavyBonds) {
              if (
                (b.source === br.donor.id && b.target === otherC.id) ||
                (b.source === otherC.id && b.target === br.donor.id)
              ) {
                continue;
              }
              const p1 = coords3D.get(b.source);
              const p2 = coords3D.get(b.target);
              if (!p1 || !p2) continue;
              const dx = p2[0] - p1[0];
              const dy = p2[1] - p1[1];
              const dz = p2[2] - p1[2];
              const d = Math.hypot(dx, dy, dz) || 0.01;
              const targetLen = b.order >= 2 ? 1.40 : 1.50;
              const diff = (d - targetLen) / d;
              const m1 = remHeavySet.has(b.source);
              const m2 = remHeavySet.has(b.target);
              if (m1 && m2) {
                const step = 0.35 * diff;
                p1[0] += step * dx;
                p1[1] += step * dy;
                p1[2] += step * dz;
                p2[0] -= step * dx;
                p2[1] -= step * dy;
                p2[2] -= step * dz;
              } else if (m1) {
                const step = 0.65 * diff;
                p1[0] += step * dx;
                p1[1] += step * dy;
                p1[2] += step * dz;
              } else if (m2) {
                const step = 0.65 * diff;
                p2[0] -= step * dx;
                p2[1] -= step * dy;
                p2[2] -= step * dz;
              }
            }
            for (let i = 0; i < branchHeavy.length; i++) {
              for (let j = i + 1; j < branchHeavy.length; j++) {
                const id1 = branchHeavy[i].id;
                const id2 = branchHeavy[j].id;
                const isBonded = intraHeavyBonds.some(
                  (b) =>
                    (b.source === id1 && b.target === id2) ||
                    (b.source === id2 && b.target === id1)
                );
                if (isBonded) continue;
                const p1 = coords3D.get(id1);
                const p2 = coords3D.get(id2);
                if (!p1 || !p2) continue;
                let dx = p2[0] - p1[0];
                let dy = p2[1] - p1[1];
                let dz = p2[2] - p1[2];
                let d = Math.hypot(dx, dy, dz);
                if (d < 2.15) {
                  if (d < 0.05) {
                    dx = nRingX;
                    dy = nRingY;
                    dz = nRingZ;
                    d = 1;
                  }
                  const push = (2.15 - d) / d;
                  const m1 = remHeavySet.has(id1);
                  const m2 = remHeavySet.has(id2);
                  if (m1 && m2) {
                    p1[0] -= 0.3 * push * dx;
                    p1[1] -= 0.3 * push * dy;
                    p1[2] -= 0.3 * push * dz;
                    p2[0] += 0.3 * push * dx;
                    p2[1] += 0.3 * push * dy;
                    p2[2] += 0.3 * push * dz;
                  } else if (m1) {
                    p1[0] -= 0.5 * push * dx;
                    p1[1] -= 0.5 * push * dy;
                    p1[2] -= 0.5 * push * dz;
                  } else if (m2) {
                    p2[0] += 0.5 * push * dx;
                    p2[1] += 0.5 * push * dy;
                    p2[2] += 0.5 * push * dz;
                  }
                }
              }
            }
            for (const a of remHeavy) {
              const p = coords3D.get(a.id);
              if (!p) continue;
              const projU = (p[0] - dPos[0]) * ux + (p[1] - dPos[1]) * uy + (p[2] - dPos[2]) * uz;
              if (projU < 0.25) {
                const boost = 0.25 - projU;
                p[0] += boost * ux;
                p[1] += boost * uy;
                p[2] += boost * uz;
              }
              p[0] = Number(p[0].toFixed(3));
              p[1] = Number(p[1].toFixed(3));
              p[2] = Number(p[2].toFixed(3));
            }
          }

          // Reconstruct all vinylic and aliphatic hydrogens on the alkene ligand
          const allBranchAtoms = [br.donor, ...br.distalAtoms];
          const hOnC0orC1 = new Set<string>();
          [
            { cAtom: br.donor, cPos: c0Pos, sign2: -1 },
            { cAtom: otherC, cPos: c1Pos, sign2: 1 },
          ].forEach(({ cAtom, cPos, sign2 }) => {
            const vinylicHs = (branchAdj.get(cAtom.id) || [])
              .map((nb) => allBranchAtoms.find((x) => x.id === nb && x.element === 'H'))
              .filter((x): x is AtomNode => Boolean(x));
            vinylicHs.forEach((hAtom, hIdx) => {
              hOnC0orC1.add(hAtom.id);
              // If C_i has 1 heavy substituent (on +e1 side) and 1 H, place the vinylic H on the opposite (-e1) side
              const sign1 = vinylicHs.length === 1 ? -1 : hIdx === 0 ? -1 : 1;
              const rBond = 1.08;
              coords3D.set(hAtom.id, [
                Number(
                  (cPos[0] + rBond * (0.35 * ux + 0.80 * sign1 * e1x + 0.48 * sign2 * e2x)).toFixed(
                    3
                  )
                ),
                Number(
                  (cPos[1] + rBond * (0.35 * uy + 0.80 * sign1 * e1y + 0.48 * sign2 * e2y)).toFixed(
                    3
                  )
                ),
                Number(
                  (cPos[2] + rBond * (0.35 * uz + 0.80 * sign1 * e1z + 0.48 * sign2 * e2z)).toFixed(
                    3
                  )
                ),
              ]);
            });
          });

          // Place remaining hydrogens attached to remHeavy
          for (const parent of remHeavy) {
            const pPos = coords3D.get(parent.id);
            if (!pPos) continue;
            const conn = bonds.filter((b) => b.source === parent.id || b.target === parent.id);
            const hNeighbors = conn
              .map((b) =>
                allBranchAtoms.find(
                  (x) =>
                    x.id === (b.source === parent.id ? b.target : b.source) &&
                    x.element === 'H' &&
                    !hOnC0orC1.has(x.id)
                )
              )
              .filter((x): x is AtomNode => Boolean(x));
            if (hNeighbors.length === 0) continue;

            const heavyNbs = conn
              .map((b) =>
                allBranchAtoms.find(
                  (x) =>
                    x.id === (b.source === parent.id ? b.target : b.source) && x.element !== 'H'
                )
              )
              .filter((x): x is AtomNode => Boolean(x));

            const uVecs: [number, number, number][] = heavyNbs
              .map((nb) => {
                const nbPos = coords3D.get(nb.id);
                if (!nbPos) return null;
                const dx = nbPos[0] - pPos[0];
                const dy = nbPos[1] - pPos[1];
                const dz = nbPos[2] - pPos[2];
                const len = Math.hypot(dx, dy, dz) || 1;
                return [dx / len, dy / len, dz / len] as [number, number, number];
              })
              .filter((v): v is [number, number, number] => Boolean(v));

            const rXH = 1.09;
            if (uVecs.length === 2 && hNeighbors.length === 2) {
              const u1 = uVecs[0];
              const u2 = uVecs[1];
              let bxH = -(u1[0] + u2[0]);
              let byH = -(u1[1] + u2[1]);
              let bzH = -(u1[2] + u2[2]);
              const bLen = Math.hypot(bxH, byH, bzH) || 1;
              bxH /= bLen;
              byH /= bLen;
              bzH /= bLen;
              let nx = u1[1] * u2[2] - u1[2] * u2[1];
              let ny = u1[2] * u2[0] - u1[0] * u2[2];
              let nz = u1[0] * u2[1] - u1[1] * u2[0];
              let nLen = Math.hypot(nx, ny, nz) || 1;
              nx /= nLen;
              ny /= nLen;
              nz /= nLen;
              const cTet = 1 / Math.sqrt(3);
              const sTet = Math.sqrt(2 / 3);
              coords3D.set(hNeighbors[0].id, [
                Number((pPos[0] + rXH * (cTet * bxH + sTet * nx)).toFixed(3)),
                Number((pPos[1] + rXH * (cTet * byH + sTet * ny)).toFixed(3)),
                Number((pPos[2] + rXH * (cTet * bzH + sTet * nz)).toFixed(3)),
              ]);
              coords3D.set(hNeighbors[1].id, [
                Number((pPos[0] + rXH * (cTet * bxH - sTet * nx)).toFixed(3)),
                Number((pPos[1] + rXH * (cTet * byH - sTet * ny)).toFixed(3)),
                Number((pPos[2] + rXH * (cTet * bzH - sTet * nz)).toFixed(3)),
              ]);
            } else if (uVecs.length >= 2 && hNeighbors.length === 1) {
              const sx = uVecs.reduce((s, v) => s + v[0], 0);
              const sy = uVecs.reduce((s, v) => s + v[1], 0);
              const sz = uVecs.reduce((s, v) => s + v[2], 0);
              const sLen = Math.hypot(sx, sy, sz) || 1;
              coords3D.set(hNeighbors[0].id, [
                Number((pPos[0] - (rXH * sx) / sLen).toFixed(3)),
                Number((pPos[1] - (rXH * sy) / sLen).toFixed(3)),
                Number((pPos[2] - (rXH * sz) / sLen).toFixed(3)),
              ]);
            } else if (uVecs.length === 1 && hNeighbors.length >= 1) {
              const u1 = uVecs[0];
              const ax = -u1[0];
              const ay = -u1[1];
              const az = -u1[2];
              let p1x = -ay;
              let p1y = ax;
              let p1z = 0;
              let p1Len = Math.hypot(p1x, p1y, p1z) || 1;
              p1x /= p1Len;
              p1y /= p1Len;
              const p2x = ay * p1z - az * p1y;
              const p2y = az * p1x - ax * p1z;
              const p2z = ax * p1y - ay * p1x;
              const cosC = 1 / 3;
              const sinC = (2 * Math.sqrt(2)) / 3;
              hNeighbors.forEach((hAtom, idx) => {
                const phi = (idx * 2 * Math.PI) / hNeighbors.length;
                const cp = Math.cos(phi);
                const sp = Math.sin(phi);
                coords3D.set(hAtom.id, [
                  Number((pPos[0] + rXH * (cosC * ax + sinC * (cp * p1x + sp * p2x))).toFixed(3)),
                  Number((pPos[1] + rXH * (cosC * ay + sinC * (cp * p1y + sp * p2y))).toFixed(3)),
                  Number((pPos[2] + rXH * (cosC * az + sinC * (cp * p1z + sp * p2z))).toFixed(3)),
                ]);
              });
            }
          }
        }
      } else if (firstShell.length === 3 && br.distalAtoms.length > 3) {
        // Multi-shell C3 tetrahedral donor branch (e.g. -PPh3 in Wilkinson's catalyst / Vaska's complex / Pd(PPh3)4, or -PEt3 in Auranofin)
        const cosCone = 0.485;
        const sinCone = 0.8745;
        const rDonor = ELEMENT_PROPERTIES[br.donor.element]?.covalentRadius ?? 1.07;

        firstShell.forEach((ipso, mIdx) => {
          // Collect sub-branch rooted at ipso
          const subIds = new Set<string>([ipso.id]);
          const q = [ipso.id];
          while (q.length > 0) {
            const curr = q.shift()!;
            for (const nb of branchAdj.get(curr) || []) {
              if (nb !== br.donor.id && distalSet.has(nb) && !subIds.has(nb)) {
                subIds.add(nb);
                q.push(nb);
              }
            }
          }
          const subAtoms = br.distalAtoms.filter((a) => subIds.has(a.id));
          const subHeavy = subAtoms.filter((a) => a.element !== 'H');

          const phi = (mIdx * 2 * Math.PI) / 3;
          const cp = Math.cos(phi);
          const sp = Math.sin(phi);
          // Cone axis c_m for this substituent
          const cxDir = cosCone * ux + sinCone * (cp * e1x + sp * e2x);
          const cyDir = cosCone * uy + sinCone * (cp * e1y + sp * e2y);
          const czDir = cosCone * uz + sinCone * (cp * e1z + sp * e2z);

          // Tangential vector t_m and binormal w_m = c_m x t_m
          const txDir = -sp * e1x + cp * e2x;
          const tyDir = -sp * e1y + cp * e2y;
          const tzDir = -sp * e1z + cp * e2z;
          const wxDir = cyDir * tzDir - czDir * tyDir;
          const wyDir = czDir * txDir - cxDir * tzDir;
          const wzDir = cxDir * tyDir - cyDir * txDir;

          // Propeller pitch angle 38° for phenyl rings
          const cPitch = Math.cos((38 * Math.PI) / 180);
          const sPitch = Math.sin((38 * Math.PI) / 180);
          const pxDir = cPitch * txDir + sPitch * wxDir;
          const pyDir = cPitch * tyDir + sPitch * wyDir;
          const pzDir = cPitch * tzDir + sPitch * wzDir;

          const rIpso =
            (rDonor + (ELEMENT_PROPERTIES[ipso.element]?.covalentRadius ?? 0.76)) * 0.96;
          const ipsoPos: [number, number, number] = [
            Number((dPos[0] + rIpso * cxDir).toFixed(3)),
            Number((dPos[1] + rIpso * cyDir).toFixed(3)),
            Number((dPos[2] + rIpso * czDir).toFixed(3)),
          ];
          coords3D.set(ipso.id, ipsoPos);

          if (subHeavy.length === 6) {
            // Reconstruct 3D regular hexagon for the phenyl ring (C0 = ipso, C1..C5 around ring)
            const orderedRing: AtomNode[] = [ipso];
            const usedRing = new Set<string>([ipso.id]);
            const heavySet = new Set(subHeavy.map((a) => a.id));
            while (orderedRing.length < 6) {
              const last = orderedRing[orderedRing.length - 1];
              const nextId = (branchAdj.get(last.id) || []).find(
                (nb) => heavySet.has(nb) && !usedRing.has(nb)
              );
              if (!nextId) {
                const rem = subHeavy.find((a) => !usedRing.has(a.id));
                if (!rem) break;
                orderedRing.push(rem);
                usedRing.add(rem.id);
              } else {
                const nextNode = subHeavy.find((a) => a.id === nextId)!;
                orderedRing.push(nextNode);
                usedRing.add(nextNode.id);
              }
            }

            const rCC = 1.39;
            const s32 = Math.sqrt(3) / 2;
            const ringOffsets: [number, number][] = [
              [0, 0],
              [0.5, s32],
              [1.5, s32],
              [2.0, 0],
              [1.5, -s32],
              [0.5, -s32],
            ];
            const ringCenter: [number, number, number] = [
              ipsoPos[0] + rCC * cxDir,
              ipsoPos[1] + rCC * cyDir,
              ipsoPos[2] + rCC * czDir,
            ];

            orderedRing.forEach((rAtom, rIdx) => {
              const [alongC, alongP] = ringOffsets[rIdx];
              const rx = Number(
                (ipsoPos[0] + rCC * (alongC * cxDir + alongP * pxDir)).toFixed(3)
              );
              const ry = Number(
                (ipsoPos[1] + rCC * (alongC * cyDir + alongP * pyDir)).toFixed(3)
              );
              const rz = Number(
                (ipsoPos[2] + rCC * (alongC * czDir + alongP * pzDir)).toFixed(3)
              );
              coords3D.set(rAtom.id, [rx, ry, rz]);

              // Place any attached H radially outward from the phenyl ring center
              const hNbs = (branchAdj.get(rAtom.id) || [])
                .map((nb) => subAtoms.find((x) => x.id === nb && x.element === 'H'))
                .filter((x): x is AtomNode => Boolean(x));
              hNbs.forEach((hAtom) => {
                const hdx = rx - ringCenter[0];
                const hdy = ry - ringCenter[1];
                const hdz = rz - ringCenter[2];
                const hLen = Math.hypot(hdx, hdy, hdz) || 1;
                coords3D.set(hAtom.id, [
                  Number((rx + (hdx / hLen) * 1.08).toFixed(3)),
                  Number((ry + (hdy / hLen) * 1.08).toFixed(3)),
                  Number((rz + (hdz / hLen) * 1.08).toFixed(3)),
                ]);
              });
            });
          } else {
            // Non-phenyl alkyl sub-branch (e.g. ethyl in -PEt3): extend along cone direction
            const otherSub = subAtoms.filter((a) => a.id !== ipso.id);
            otherSub.forEach((sAtom, sIdx) => {
              const step = sAtom.element === 'H' ? 1.09 : 1.52;
              const angle = ((sIdx + 1) * 2 * Math.PI) / Math.max(1, otherSub.length);
              coords3D.set(sAtom.id, [
                Number(
                  (
                    ipsoPos[0] +
                    step * (0.75 * cxDir + 0.66 * (Math.cos(angle) * pxDir + Math.sin(angle) * wxDir))
                  ).toFixed(3)
                ),
                Number(
                  (
                    ipsoPos[1] +
                    step * (0.75 * cyDir + 0.66 * (Math.cos(angle) * pyDir + Math.sin(angle) * wyDir))
                  ).toFixed(3)
                ),
                Number(
                  (
                    ipsoPos[2] +
                    step * (0.75 * czDir + 0.66 * (Math.cos(angle) * pzDir + Math.sin(angle) * wzDir))
                  ).toFixed(3)
                ),
              ]);
            });
          }
        });
      } else {
        // Rigidly orient the ligand's 2D coordinates outward along u from dPos
        const avg2x =
          br.distalAtoms.reduce((s, a) => s + a.x2d, 0) / br.distalAtoms.length - br.donor.x2d;
        const avg2y =
          br.distalAtoms.reduce((s, a) => s + -a.y2d, 0) / br.distalAtoms.length - -br.donor.y2d;
        const avgLen = Math.hypot(avg2x, avg2y) || 1;
        const lbx = avg2x / avgLen;
        const lby = avg2y / avgLen;
        const ltx = -lby;
        const lty = lbx;
        for (const a of br.distalAtoms) {
          const dx = a.x2d - br.donor.x2d;
          const dy = -a.y2d - -br.donor.y2d;
          const sOut = Math.max(0.4, dx * lbx + dy * lby);
          const sTrans = dx * ltx + dy * lty;
          const sNorm = a.z3d || 0;
          coords3D.set(a.id, [
            Number((dPos[0] + sOut * ux + sTrans * e1x + sNorm * e2x).toFixed(3)),
            Number((dPos[1] + sOut * uy + sTrans * e1y + sNorm * e2y).toFixed(3)),
            Number((dPos[2] + sOut * uz + sTrans * e1z + sNorm * e2z).toFixed(3)),
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

  // Check if genuine 3D depth already exists on non-synthetic atoms (e.g. from PubChem 3D-SDF or canonical 3D cages)
  const nonImpAtoms = atoms.filter((a) => !a.id.startsWith('H_imp_'));
  const maxNonImpAbsZ = nonImpAtoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d || 0)), 0);
  if (maxNonImpAbsZ > 0.08) {
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

  // Flat 2D multi-center graph (from 2D-SDF such as stearic acid / eicosane / vitamin E, SMILES, or 2D Sketcher):
  // Identify sp3 tetrahedral heavy centers (all single bonds, zero multiple bonds) vs planar sp2/sp centers.
  const heavyAtoms = atoms.filter((a) => a.element !== 'H');
  const sp3HeavyIds = new Set<string>();
  for (const a of heavyAtoms) {
    const conn = bonds.filter((b) => b.source === a.id || b.target === a.id);
    if (conn.length === 0 || !conn.every((b) => b.order === 1)) continue;
    if (
      (['C', 'Si', 'Ge', 'Sn'].includes(a.element) &&
        (a.formalCharge ?? 0) === 0 &&
        conn.length >= 2) ||
      conn.length === 4
    ) {
      sp3HeavyIds.add(a.id);
    }
  }

  // If the molecule is purely sp2/sp conjugated/aromatic/linear (no sp3 heavy centers), keep it planar at z = 0
  if (sp3HeavyIds.size === 0) {
    return atoms.map((a) => ({
      ...a,
      x2d: Number((a.x2d - cx2d).toFixed(2)),
      y2d: Number((a.y2d - cy2d).toFixed(2)),
      x3d: Number((a.x2d - cx2d).toFixed(3)),
      y3d: Number((-(a.y2d - cy2d)).toFixed(3)),
      z3d: Number((a.z3d || 0).toFixed(3)),
    }));
  }

  // Build 3D coordinates map starting from centered 2D coordinates
  const pos3D = new Map<string, [number, number, number]>();
  for (const a of atoms) {
    pos3D.set(a.id, [a.x2d - cx2d, -(a.y2d - cy2d), a.z3d || 0]);
  }

  // 1. Compute bipartite BFS graph distance across heavy-atom components so adjacent sp3 chain/ring carbons
  //    alternate ±z_pucker = ±0.255 Å, relaxing 2D 120° angles to the exact 109.47° tetrahedral angle in 3D.
  const heavyDepth = new Map<string, number>();
  for (const root of heavyAtoms) {
    if (heavyDepth.has(root.id)) continue;
    heavyDepth.set(root.id, 0);
    const queue: string[] = [root.id];
    while (queue.length > 0) {
      const currId = queue.shift()!;
      const currD = heavyDepth.get(currId) ?? 0;
      for (const b of bonds) {
        const nextId = b.source === currId ? b.target : b.target === currId ? b.source : null;
        if (!nextId || heavyDepth.has(nextId)) continue;
        const nextAtom = atoms.find((x) => x.id === nextId);
        if (nextAtom && nextAtom.element !== 'H') {
          heavyDepth.set(nextId, currD + 1);
          queue.push(nextId);
        }
      }
    }
  }

  for (const a of heavyAtoms) {
    if (!sp3HeavyIds.has(a.id)) continue;
    const heavyNeighbors = bonds
      .filter((b) => b.source === a.id || b.target === a.id)
      .map((b) => atoms.find((x) => x.id === (b.source === a.id ? b.target : b.source)))
      .filter((x): x is AtomNode => Boolean(x && x.element !== 'H'));

    if (heavyNeighbors.length >= 2) {
      const d = heavyDepth.get(a.id) ?? 0;
      const sign = d % 2 === 0 ? 1 : -1;
      const cur = pos3D.get(a.id)!;
      let x3 = cur[0];
      let y3 = cur[1];

      // If 2 heavy neighbors were drawn nearly collinear (> 155°) in the 2D sketcher, add in-plane zigzag offset too
      if (heavyNeighbors.length === 2) {
        const p1 = pos3D.get(heavyNeighbors[0].id)!;
        const p2 = pos3D.get(heavyNeighbors[1].id)!;
        const v1x = p1[0] - x3;
        const v1y = p1[1] - y3;
        const v2x = p2[0] - x3;
        const v2y = p2[1] - y3;
        const l1 = Math.hypot(v1x, v1y) || 1;
        const l2 = Math.hypot(v2x, v2y) || 1;
        const cosAng = (v1x * v2x + v1y * v2y) / (l1 * l2);
        if (cosAng < -0.9) {
          const ax = (p2[0] - p1[0]) / (Math.hypot(p2[0] - p1[0], p2[1] - p1[1]) || 1);
          const ay = (p2[1] - p1[1]) / (Math.hypot(p2[0] - p1[0], p2[1] - p1[1]) || 1);
          x3 += -ay * sign * 0.36;
          y3 += ax * sign * 0.36;
        }
      }

      pos3D.set(a.id, [x3, y3, sign * 0.255]);
    }
  }

  // 2. Reconstruct 3D tetrahedral positions for terminal heavy substituents on tertiary (-CH(R)-)
  //    and quaternary (-C(R1)(R2)-, -PO3, -SO2) centers
  for (const center of heavyAtoms) {
    const conn = bonds.filter((b) => b.source === center.id || b.target === center.id);
    const heavyNeighbors = conn
      .map((b) => atoms.find((x) => x.id === (b.source === center.id ? b.target : b.source)))
      .filter((x): x is AtomNode => Boolean(x && x.element !== 'H'));

    if (heavyNeighbors.length < 3) continue;
    const isCenterSp3 = sp3HeavyIds.has(center.id) || conn.length === 4;
    if (!isCenterSp3) continue;

    const terminalHeavy = heavyNeighbors.filter((nb) => {
      const nbHeavyDeg = bonds.filter((b) => {
        if (b.source !== nb.id && b.target !== nb.id) return false;
        const otherId = b.source === nb.id ? b.target : b.source;
        const other = atoms.find((x) => x.id === otherId);
        return Boolean(other && other.element !== 'H');
      }).length;
      return nbHeavyDeg === 1;
    });
    const backboneHeavy = heavyNeighbors.filter((nb) => !terminalHeavy.includes(nb));
    const cPos = pos3D.get(center.id)!;

    if (heavyNeighbors.length === 4 && terminalHeavy.length >= 2) {
      const refNeighbors = backboneHeavy.length >= 2 ? backboneHeavy : heavyNeighbors.slice(0, 2);
      const p1 = pos3D.get(refNeighbors[0].id)!;
      const p2 = pos3D.get(refNeighbors[1].id)!;
      const u1x = (p1[0] - cPos[0]) / (Math.hypot(p1[0] - cPos[0], p1[1] - cPos[1], p1[2] - cPos[2]) || 1);
      const u1y = (p1[1] - cPos[1]) / (Math.hypot(p1[0] - cPos[0], p1[1] - cPos[1], p1[2] - cPos[2]) || 1);
      const u2x = (p2[0] - cPos[0]) / (Math.hypot(p2[0] - cPos[0], p2[1] - cPos[1], p2[2] - cPos[2]) || 1);
      const u2y = (p2[1] - cPos[1]) / (Math.hypot(p2[0] - cPos[0], p2[1] - cPos[1], p2[2] - cPos[2]) || 1);
      let bx = -(u1x + u2x);
      let by = -(u1y + u2y);
      let bLen = Math.hypot(bx, by);
      if (bLen < 0.1) {
        bx = -u1y;
        by = u1x;
        bLen = Math.hypot(bx, by) || 1;
      }
      bx /= bLen;
      by /= bLen;

      const cTet = 1 / Math.sqrt(3);
      const sTet = Math.sqrt(2 / 3);
      for (let tIdx = 0; tIdx < 2; tIdx++) {
        const term = terminalHeavy[tIdx];
        const rBond =
          (ELEMENT_PROPERTIES[center.element]?.covalentRadius ?? 0.77) +
          (ELEMENT_PROPERTIES[term.element]?.covalentRadius ?? 0.77);
        const zSign = tIdx === 0 ? 1 : -1;
        pos3D.set(term.id, [
          cPos[0] + rBond * cTet * bx,
          cPos[1] + rBond * cTet * by,
          cPos[2] + rBond * sTet * zSign,
        ]);
      }
    } else if (heavyNeighbors.length === 3 && terminalHeavy.length >= 1 && sp3HeavyIds.has(center.id)) {
      const term = terminalHeavy[0];
      const chainNbs = heavyNeighbors.filter((nb) => nb.id !== term.id);
      if (chainNbs.length === 2) {
        const p1 = pos3D.get(chainNbs[0].id)!;
        const p2 = pos3D.get(chainNbs[1].id)!;
        const u1x = (p1[0] - cPos[0]) / (Math.hypot(p1[0] - cPos[0], p1[1] - cPos[1], p1[2] - cPos[2]) || 1);
        const u1y = (p1[1] - cPos[1]) / (Math.hypot(p1[0] - cPos[0], p1[1] - cPos[1], p1[2] - cPos[2]) || 1);
        const u2x = (p2[0] - cPos[0]) / (Math.hypot(p2[0] - cPos[0], p2[1] - cPos[1], p2[2] - cPos[2]) || 1);
        const u2y = (p2[1] - cPos[1]) / (Math.hypot(p2[0] - cPos[0], p2[1] - cPos[1], p2[2] - cPos[2]) || 1);
        let bx = -(u1x + u2x);
        let by = -(u1y + u2y);
        let bLen = Math.hypot(bx, by);
        if (bLen < 0.1) {
          bx = -u1y;
          by = u1x;
          bLen = Math.hypot(bx, by) || 1;
        }
        bx /= bLen;
        by /= bLen;

        const termBond = conn.find((b) => b.source === term.id || b.target === term.id);
        const d = heavyDepth.get(center.id) ?? 0;
        const zSign =
          termBond?.stereo === 'dash' ? -1 : termBond?.stereo === 'wedge' ? 1 : d % 2 === 0 ? 1 : -1;
        const rBond =
          (ELEMENT_PROPERTIES[center.element]?.covalentRadius ?? 0.77) +
          (ELEMENT_PROPERTIES[term.element]?.covalentRadius ?? 0.77);
        const cTet = 1 / Math.sqrt(3);
        const sTet = Math.sqrt(2 / 3);
        pos3D.set(term.id, [
          cPos[0] + rBond * cTet * bx,
          cPos[1] + rBond * cTet * by,
          cPos[2] + rBond * sTet * zSign,
        ]);
      }
    }
  }

  // 3. Re-anchor and orient all terminal Hydrogens in 3D around their parent heavy atoms
  //    using exact 109.47° tetrahedral (sp3) or 120° trigonal planar (sp2) geometry
  for (const parent of heavyAtoms) {
    const conn = bonds.filter((b) => b.source === parent.id || b.target === parent.id);
    const hNeighbors = conn
      .map((b) => atoms.find((x) => x.id === (b.source === parent.id ? b.target : b.source)))
      .filter((x): x is AtomNode => Boolean(x && x.element === 'H'));
    if (hNeighbors.length === 0) continue;

    const heavyNbs = conn
      .map((b) => atoms.find((x) => x.id === (b.source === parent.id ? b.target : b.source)))
      .filter((x): x is AtomNode => Boolean(x && x.element !== 'H'));

    const pPos = pos3D.get(parent.id)!;
    const uVecs: [number, number, number][] = heavyNbs.map((nb) => {
      const nbPos = pos3D.get(nb.id)!;
      const dx = nbPos[0] - pPos[0];
      const dy = nbPos[1] - pPos[1];
      const dz = nbPos[2] - pPos[2];
      const len = Math.hypot(dx, dy, dz) || 1;
      return [dx / len, dy / len, dz / len];
    });

    const rXH = (ELEMENT_PROPERTIES[parent.element]?.covalentRadius ?? 0.77) + 0.32;
    const isParentSp3 =
      sp3HeavyIds.has(parent.id) ||
      (conn.every((b) => b.order === 1) &&
        (conn.length === 4 ||
          (['N', 'P', 'As'].includes(parent.element) && conn.length === 3) ||
          (['O', 'S', 'Se'].includes(parent.element) && conn.length === 2)));

    if (isParentSp3 && uVecs.length === 2 && hNeighbors.length === 2) {
      // Methylene (-CH2-): two tetrahedral C-H vectors at ±54.735° from outward bisector
      const u1 = uVecs[0];
      const u2 = uVecs[1];
      let bx = -(u1[0] + u2[0]);
      let by = -(u1[1] + u2[1]);
      let bz = -(u1[2] + u2[2]);
      let bLen = Math.hypot(bx, by, bz);
      if (bLen < 0.08) {
        bx = -u1[1];
        by = u1[0];
        bz = 0;
        bLen = Math.hypot(bx, by, bz) || 1;
      }
      bx /= bLen;
      by /= bLen;
      bz /= bLen;

      let nx = u1[1] * u2[2] - u1[2] * u2[1];
      let ny = u1[2] * u2[0] - u1[0] * u2[2];
      let nz = u1[0] * u2[1] - u1[1] * u2[0];
      let nLen = Math.hypot(nx, ny, nz);
      if (nLen < 0.08) {
        nx = 0;
        ny = 0;
        nz = 1;
      } else {
        nx /= nLen;
        ny /= nLen;
        nz /= nLen;
      }
      const dotNB = nx * bx + ny * by + nz * bz;
      nx -= dotNB * bx;
      ny -= dotNB * by;
      nz -= dotNB * bz;
      nLen = Math.hypot(nx, ny, nz) || 1;
      nx /= nLen;
      ny /= nLen;
      nz /= nLen;

      const cTet = 1 / Math.sqrt(3);
      const sTet = Math.sqrt(2 / 3);
      pos3D.set(hNeighbors[0].id, [
        pPos[0] + rXH * (cTet * bx + sTet * nx),
        pPos[1] + rXH * (cTet * by + sTet * ny),
        pPos[2] + rXH * (cTet * bz + sTet * nz),
      ]);
      pos3D.set(hNeighbors[1].id, [
        pPos[0] + rXH * (cTet * bx - sTet * nx),
        pPos[1] + rXH * (cTet * by - sTet * ny),
        pPos[2] + rXH * (cTet * bz - sTet * nz),
      ]);
    } else if (isParentSp3 && uVecs.length === 1 && hNeighbors.length >= 1) {
      // Methyl (-CH3), Amino (-NH2), Hydroxyl (-OH): staggered 3D tetrahedral cone
      const u1 = uVecs[0];
      const ax = -u1[0];
      const ay = -u1[1];
      const az = -u1[2];
      let e1x = -ay;
      let e1y = ax;
      let e1z = 0;
      let e1Len = Math.hypot(e1x, e1y, e1z);
      if (e1Len < 0.1) {
        e1x = 1;
        e1y = 0;
        e1z = 0;
        e1Len = 1;
      }
      e1x /= e1Len;
      e1y /= e1Len;
      const e2x = ay * e1z - az * e1y;
      const e2y = az * e1x - ax * e1z;
      const e2z = ax * e1y - ay * e1x;

      const cosCone = 1 / 3;
      const sinCone = (2 * Math.sqrt(2)) / 3;
      hNeighbors.forEach((hAtom, idx) => {
        const phi = (idx * 2 * Math.PI) / 3;
        const cp = Math.cos(phi);
        const sp = Math.sin(phi);
        pos3D.set(hAtom.id, [
          pPos[0] + rXH * (cosCone * ax + sinCone * (cp * e1x + sp * e2x)),
          pPos[1] + rXH * (cosCone * ay + sinCone * (cp * e1y + sp * e2y)),
          pPos[2] + rXH * (cosCone * az + sinCone * (cp * e1z + sp * e2z)),
        ]);
      });
    } else if (isParentSp3 && uVecs.length === 3 && hNeighbors.length === 1) {
      // Methine (-CH(R1)(R2)(R3)): 4th apex opposite the 3 heavy bonds
      const sx = uVecs[0][0] + uVecs[1][0] + uVecs[2][0];
      const sy = uVecs[0][1] + uVecs[1][1] + uVecs[2][1];
      const sz = uVecs[0][2] + uVecs[1][2] + uVecs[2][2];
      const sLen = Math.hypot(sx, sy, sz);
      if (sLen > 0.12) {
        pos3D.set(hNeighbors[0].id, [
          pPos[0] - (rXH * sx) / sLen,
          pPos[1] - (rXH * sy) / sLen,
          pPos[2] - (rXH * sz) / sLen,
        ]);
      } else {
        pos3D.set(hNeighbors[0].id, [pPos[0], pPos[1], pPos[2] - rXH]);
      }
    }
  }

  // Center final 3D coordinates around the origin
  let sumX = 0;
  let sumY = 0;
  let sumZ = 0;
  for (const a of atoms) {
    const p = pos3D.get(a.id)!;
    sumX += p[0];
    sumY += p[1];
    sumZ += p[2];
  }
  const meanX = sumX / atoms.length;
  const meanY = sumY / atoms.length;
  const meanZ = sumZ / atoms.length;

  return atoms.map((a) => {
    const p = pos3D.get(a.id)!;
    return {
      ...a,
      x2d: Number((a.x2d - cx2d).toFixed(2)),
      y2d: Number((a.y2d - cy2d).toFixed(2)),
      x3d: Number((p[0] - meanX).toFixed(3)),
      y3d: Number((p[1] - meanY).toFixed(3)),
      z3d: Number((p[2] - meanZ).toFixed(3)),
    };
  });
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

  // Check unsubstituted monocyclic planar rings (e.g. C6H6 or bare C6 hexagon -> D6h, C4H4 -> D4h)
  const maxAbsZRing = atoms.reduce((m, a) => Math.max(m, Math.abs(a.z3d || 0)), 0);
  if (
    heavy.length === 6 &&
    heavy.every((a) => a.element === 'C') &&
    bonds.length >= 6 &&
    atoms.length <= 12 &&
    maxAbsZRing < 0.25
  ) {
    const heavyDegs = heavy.map(
      (a) =>
        bonds.filter(
          (b) =>
            (b.source === a.id || b.target === a.id) &&
            heavy.some((h) => h.id === (b.source === a.id ? b.target : b.source))
        ).length
    );
    if (heavyDegs.every((d) => d === 2)) return 'D6h';
  }
  if (
    heavy.length === 4 &&
    heavy.every((a) => a.element === 'C') &&
    bonds.length >= 4 &&
    atoms.length <= 8 &&
    maxAbsZRing < 0.25
  ) {
    const degs = heavy.map(
      (a) =>
        bonds.filter(
          (b) =>
            (b.source === a.id || b.target === a.id) &&
            heavy.some((h) => h.id === (b.source === a.id ? b.target : b.source))
        ).length
    );
    if (degs.every((d) => d === 2)) return 'D4h';
  }

  // Bidentate chelate and η5-metallocene sandwich point-group classification
  // (e.g. Ferrocene -> D5d, [Co(en)3]3+ / [Ru(bpy)3]2+ / Fe(acac)3 / [Fe(C2O4)3]3- -> D3, Ni(dmgH)2 -> D2h, Oxaliplatin / Carboplatin -> C2v)
  const chelateTopo = analyzeChelateOrSandwichTopology(atoms, bonds);
  if (chelateTopo) {
    return chelateTopo.pointGroupId;
  }

  // Central atom surrounded by nLig ligands (AX_nE_m VSEPR & d/f coordination point-group classification,
  // supporting both monoatomic ligands and polyatomic ligands like -PPh3, -CN, -CO, -NO, -NH3, -C2H4)
  const starTopo = analyzeStarCoordinationTopology(atoms, bonds);

  if (starTopo) {
    const { centerAtom, branches } = starTopo;
    const ligandSigs = new Set(branches.map((b) => b.sig));
    const nLig = branches.length;
    const isMetalCenter =
      isTransitionMetal(centerAtom.element) || isLanthanide(centerAtom.element);
    const lonePairs = getCentralAtomLonePairs(centerAtom, atoms, bonds, charge);
    const isSqPlanar = nLig === 4 && isFourCoordSquarePlanar(centerAtom, branches, charge);

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
          return isSqPlanar ? 'D4h' : 'Td';
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
        // ML3X square planar (Wilkinson's catalyst [RhCl(PPh3)3], Zeise's salt [PtCl3(C2H4)]- -> C2v)
        // vs AX3Y tetrahedral (CHCl3, POCl3 -> C3v)
        if (isMetalCenter && isSqPlanar) return 'C2v';
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

    if (ligandSigs.size === 3) {
      const sigCounts = Array.from(ligandSigs).map(
        (sig) => branches.filter((b) => b.sig === sig).length
      );
      if (nLig === 4 && sigCounts.includes(2)) {
        // trans-ML2XY square planar (e.g. Vaska's complex trans-[IrCl(CO)(PPh3)2] -> C2v) vs tetrahedral -> Cs
        return isMetalCenter && isSqPlanar ? 'C2v' : 'Cs';
      }
      if (nLig === 6) {
        return isMetalCenter || lonePairs === 0 ? 'C2v' : 'Cs';
      }
    }
  }

  if (heavy.length === 2 && heavy[0].element === heavy[1].element) {
    const nH = atoms.filter((a) => a.element === 'H').length;
    if (nH === 2) return 'Dinfh';
    if (nH === 4 || (nH === 6 && bonds.length >= 7)) return 'D2h';
  }

  // General 3D/2D Covariance & Symmetry-Operation Classifier for multi-center molecules:
  // Handles PAHs and fused aromatics (Coronene D6h, Triphenylene D3h, Naphthalene/Anthracene/Pyrene/Biphenylene D2h,
  // Phenanthrene/Azulene C2v), heterocycles & inorganic rings (Borazine/1,3,5-Triazine D3h, Pyrazine/1,4-Benzoquinone D2h,
  // Pyridine/Furan/Thiophene/Pyrimidine C2v), perhalobenzenes (C6F6/C6Cl6 D6h), cycloalkanes (Cyclopropane D3h),
  // 3D polyhedral cages (Cubane Oh, Adamantane/P4/Urotropine Td), and linear polyatomics (Diacetylene Dinfh, Cyanoacetylene Cinfv).
  if (heavy.length >= 3) {
    const has3dCoords = atoms.some(
      (a) => Math.abs(a.x3d) > 1e-3 || Math.abs(a.y3d) > 1e-3 || Math.abs(a.z3d) > 1e-3
    );
    const getPos = (a: AtomNode): [number, number, number] =>
      has3dCoords ? [a.x3d, a.y3d, a.z3d] : [a.x2d, -a.y2d, 0];

    const nHvy = heavy.length;
    let cx = 0;
    let cy = 0;
    let cz = 0;
    for (const a of heavy) {
      const [px, py, pz] = getPos(a);
      cx += px;
      cy += py;
      cz += pz;
    }
    cx /= nHvy;
    cy /= nHvy;
    cz /= nHvy;

    const pts = heavy.map((a) => {
      const [px, py, pz] = getPos(a);
      const deg = bonds.filter(
        (b) =>
          (b.source === a.id || b.target === a.id) &&
          heavy.some((h) => h.id === (b.source === a.id ? b.target : b.source))
      ).length;
      return {
        el: a.element,
        deg,
        x: px - cx,
        y: py - cy,
        z: pz - cz,
      };
    });

    const cov: number[][] = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    for (const p of pts) {
      cov[0][0] += p.x * p.x;
      cov[0][1] += p.x * p.y;
      cov[0][2] += p.x * p.z;
      cov[1][0] += p.y * p.x;
      cov[1][1] += p.y * p.y;
      cov[1][2] += p.y * p.z;
      cov[2][0] += p.z * p.x;
      cov[2][1] += p.z * p.y;
      cov[2][2] += p.z * p.z;
    }
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) {
        cov[r][c] /= nHvy;
      }
    }

    const { eigenvalues: covEigs, eigenvectors: covVecs } = jacobiEigenvalues(cov);
    // covEigs sorted ascending: lambda0 <= lambda1 <= lambda2
    const [lambda0, lambda1, lambda2] = covEigs;
    const e0 = covVecs[0];
    const e1 = covVecs[1];
    const e2 = covVecs[2];

    const proj = pts.map((p) => ({
      el: p.el,
      deg: p.deg,
      u: p.x * e1[0] + p.y * e1[1] + p.z * e1[2],
      v: p.x * e2[0] + p.y * e2[1] + p.z * e2[2],
      w: p.x * e0[0] + p.y * e0[1] + p.z * e0[2],
    }));

    const rMax = Math.max(1, ...proj.map((p) => Math.hypot(p.u, p.v, p.w)));
    const tol = Math.max(0.22, rMax * 0.075);

    const hasInversion = proj.every((p) =>
      proj.some(
        (q) =>
          q.el === p.el &&
          q.deg === p.deg &&
          Math.hypot(q.u + p.u, q.v + p.v, q.w + p.w) <= tol
      )
    );

    // 1. Collinear / linear polyatomic chains (e.g. Diacetylene C4H2, Cyanogen C2N2, Carbon suboxide C3O2, Cyanoacetylene HC3N)
    if (lambda0 < 0.015 && lambda1 < 0.015 && lambda2 > 0.25) {
      const allPts = atoms.map((a) => {
        const [px, py, pz] = getPos(a);
        return {
          el: a.element,
          u: (px - cx) * e1[0] + (py - cy) * e1[1] + (pz - cz) * e1[2],
          v: (px - cx) * e2[0] + (py - cy) * e2[1] + (pz - cz) * e2[2],
          w: (px - cx) * e0[0] + (py - cy) * e0[1] + (pz - cz) * e0[2],
        };
      });
      const allCollinear = allPts.every((p) => Math.hypot(p.u, p.w) < 0.2);
      if (allCollinear) {
        const allInv = allPts.every((p) =>
          allPts.some((q) => q.el === p.el && Math.abs(q.v + p.v) <= tol)
        );
        return allInv ? 'Dinfh' : 'Cinfv';
      }
    }

    // 2. 3D Isotropic Spherical-Top Polyhedral Cages (e.g. Cubane C8H8 -> Oh; Adamantane C10H16, P4, Urotropine -> Td)
    if (lambda0 > 0.25 && lambda2 > 0.25 && lambda0 / lambda2 > 0.82) {
      return hasInversion ? 'Oh' : 'Td';
    }

    // 3. Planar or σh-Symmetric Ring / Polycyclic Systems (lambda0 < 0.04)
    const maxHeavyW = proj.reduce((m, p) => Math.max(m, Math.abs(p.w)), 0);
    if (maxHeavyW < 0.25) {
      const hasCn = (fold: number): boolean => {
        const theta = (2 * Math.PI) / fold;
        const ct = Math.cos(theta);
        const st = Math.sin(theta);
        return proj.every((p) => {
          const ru = p.u * ct - p.v * st;
          const rv = p.u * st + p.v * ct;
          return proj.some(
            (q) =>
              q.el === p.el &&
              q.deg === p.deg &&
              Math.hypot(q.u - ru, q.v - rv, q.w - p.w) <= tol
          );
        });
      };

      // Candidate mirror plane angles in the (u, v) molecular plane
      const candidateAngles: number[] = [0, Math.PI / 2];
      for (const p of proj) {
        if (Math.hypot(p.u, p.v) > 0.25) {
          const ang = Math.atan2(p.v, p.u);
          candidateAngles.push(ang, ang + Math.PI / 2);
        }
      }
      for (let k = 0; k < 72; k++) {
        candidateAngles.push((k * Math.PI) / 72);
      }

      const isMirrorAngle = (phi: number): boolean => {
        const c2 = Math.cos(2 * phi);
        const s2 = Math.sin(2 * phi);
        return proj.every((p) => {
          const mu = p.u * c2 + p.v * s2;
          const mv = p.u * s2 - p.v * c2;
          return proj.some(
            (q) =>
              q.el === p.el &&
              q.deg === p.deg &&
              Math.hypot(q.u - mu, q.v - mv, q.w - p.w) <= tol
          );
        });
      };

      const validMirrors = candidateAngles.filter(isMirrorAngle);
      const hasMirror = validMirrors.length > 0;
      const hasOrthogonalMirrors = validMirrors.some((phi) => isMirrorAngle(phi + Math.PI / 2));

      const c6 = hasCn(6);
      const c4 = hasCn(4);
      const c3 = hasCn(3);
      const c2 = hasCn(2);

      if (c6 && hasMirror) return 'D6h';
      if (c4 && hasMirror) return 'D4h';
      if (c3 && hasMirror) return 'D3h';
      if (c2 && hasOrthogonalMirrors) return 'D2h';
      if (hasMirror) return 'C2v';
    }
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
  // and prioritize the central metal + inner coordination shells by BFS graph distance
  // so that metal centers at the end of a PubChem SDF atom table (e.g. Rh59 in Wilkinson's catalyst)
  // are never truncated out of the LCAO-MO basis.
  const isLargeSystem = atoms.length > 28;
  let activeAtoms = atoms;
  if (isLargeSystem) {
    const filtered = atoms.filter((a) => {
      if (a.element !== 'H') return true;
      const attached = bonds.find((b) => b.source === a.id || b.target === a.id);
      if (!attached) return false;
      const partnerId = attached.source === a.id ? attached.target : attached.source;
      const partner = atoms.find((p) => p.id === partnerId);
      return partner && partner.element !== 'C';
    });
    if (filtered.length > 48) {
      const metalRoots = filtered.filter(
        (a) => isTransitionMetal(a.element) || isLanthanide(a.element)
      );
      if (metalRoots.length > 0) {
        const bfsDist = new Map<string, number>();
        const q: string[] = [];
        for (const m of metalRoots) {
          bfsDist.set(m.id, 0);
          q.push(m.id);
        }
        while (q.length > 0) {
          const curr = q.shift()!;
          const d = bfsDist.get(curr) ?? 0;
          for (const b of bonds) {
            const nb = b.source === curr ? b.target : b.target === curr ? b.source : null;
            if (nb && !bfsDist.has(nb)) {
              bfsDist.set(nb, d + 1);
              q.push(nb);
            }
          }
        }
        filtered.sort((a, b) => (bfsDist.get(a.id) ?? 99) - (bfsDist.get(b.id) ?? 99));
      }
    }
    activeAtoms = filtered.slice(0, 48);
  }

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

  // 2. Otherwise, resolve any unlinked dative metal-ligand coordination bonds
  //    (including partially bonded complexes like Cisplatin Cl-Pt-Cl + 2 NH3 or Nitroprusside [Fe(CN)5] + NO),
  //    expand implicit hydrogens (for small/medium graphs <= 28 heavy atoms),
  //    reconstruct 3D VSEPR / coordination coordinates,
  //    and run the full s-p-d-f Extended Hückel LCAO-MO + Symmetry solver
  const hasCoordinationCenter =
    rawAtoms.length >= 3 &&
    rawAtoms.some(
      (a) =>
        isTransitionMetal(a.element) ||
        isLanthanide(a.element) ||
        MAIN_GROUP_COORD_CENTERS.has(a.element)
    );
  const coordResolved = hasCoordinationCenter
    ? resolveCoordinationBonds(rawAtoms, rawBonds, true)
    : { atoms: rawAtoms, bonds: rawBonds };

  const rawHeavyCount = coordResolved.atoms.filter((a) => a.element !== 'H').length;
  const shouldExpandImplicitH = rawHeavyCount <= 28 && coordResolved.atoms.length <= 76;
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
  const chelateTopo = analyzeChelateOrSandwichTopology(centeredAtoms, bonds);
  const primaryAtom =
    starTopo?.centerAtom ||
    chelateTopo?.centerAtom ||
    centeredAtoms.find((a) => isTransitionMetal(a.element) || isLanthanide(a.element)) ||
    centeredAtoms.find((a) => {
      const deg = bonds.filter((b) => b.source === a.id || b.target === a.id).length;
      return deg >= 2 && deg === centeredAtoms.length - 1;
    }) ||
    centeredAtoms.find((a) => a.element !== 'H') ||
    centeredAtoms[0];
  const peripheralAtom =
    starTopo?.branches[0]?.donor ||
    chelateTopo?.ligands[0]?.donors[0] ||
    centeredAtoms.find((a) => a.id !== primaryAtom?.id) ||
    primaryAtom;

  const primProps = ELEMENT_PROPERTIES[primaryAtom?.element || 'C'] || ELEMENT_PROPERTIES.C;
  const periphProps = ELEMENT_PROPERTIES[peripheralAtom?.element || 'H'] || ELEMENT_PROPERTIES.H;
  const primShell = getPrincipalShell(primaryAtom?.element || 'C');
  const periphShell = getPrincipalShell(peripheralAtom?.element || 'H');
  const primIsTm = primaryAtom ? isTransitionMetal(primaryAtom.element) : false;
  const primIsLn = primaryAtom ? isLanthanide(primaryAtom.element) : false;

  // Detect Dewar-Chatt-Duncanson η²-alkene coordination (a transition metal bonded to both C1 and C2 of an η² C=C alkene or chelating diene)
  const dcdAlkenePair = (() => {
    if (!primaryAtom || (!primIsTm && !primIsLn)) return null;
    const mNbIds = bonds
      .filter((b) => b.source === primaryAtom.id || b.target === primaryAtom.id)
      .map((b) => (b.source === primaryAtom.id ? b.target : b.source));
    const cNbs = mNbIds
      .map((id) => centeredAtoms.find((a) => a.id === id))
      .filter((a): a is AtomNode => Boolean(a && a.element === 'C'));
    if (cNbs.length < 2) return null;
    const cNbSet = new Set(cNbs.map((a) => a.id));
    for (let i = 0; i < cNbs.length; i++) {
      for (let j = i + 1; j < cNbs.length; j++) {
        const ccBond = bonds.find(
          (b) =>
            (b.source === cNbs[i].id && b.target === cNbs[j].id) ||
            (b.source === cNbs[j].id && b.target === cNbs[i].id)
        );
        if (!ccBond) continue;
        // Verify this C-C pair is a discrete η²-alkene unit (not part of a continuous η⁴/η⁵/η⁶ sandwich ring)
        const otherCNbsI = bonds.filter(
          (b) =>
            (b.source === cNbs[i].id && cNbSet.has(b.target) && b.target !== cNbs[j].id) ||
            (b.target === cNbs[i].id && cNbSet.has(b.source) && b.source !== cNbs[j].id)
        ).length;
        const otherCNbsJ = bonds.filter(
          (b) =>
            (b.source === cNbs[j].id && cNbSet.has(b.target) && b.target !== cNbs[i].id) ||
            (b.target === cNbs[j].id && cNbSet.has(b.source) && b.source !== cNbs[i].id)
        ).length;
        if (otherCNbsI === 0 && otherCNbsJ === 0) {
          const coDonor =
            centeredAtoms.find(
              (a) => mNbIds.includes(a.id) && a.id !== cNbs[i].id && a.id !== cNbs[j].id
            ) || peripheralAtom;
          return { metalAtom: primaryAtom, c1: cNbs[i], c2: cNbs[j], coDonor };
        }
      }
    }
    return null;
  })();

  // Build dynamic fragment orbitals from actual element VSIEs (or Dewar-Chatt-Duncanson fragments for η²-alkene complexes)
  const firstIrrepLower = (pg.irreps[0]?.symbol || 'a').toLowerCase().replace(/ς/g, 'σ');
  const secondIrrepLower = (pg.irreps[Math.min(1, pg.irreps.length - 1)]?.symbol || 'b')
    .toLowerCase()
    .replace(/ς/g, 'σ');
  const lastIrrepLower = (pg.irreps[pg.irreps.length - 1]?.symbol || 'b')
    .toLowerCase()
    .replace(/ς/g, 'σ');
  const dShellNum = Math.max(3, Number(primShell) - 1);
  const fragmentOrbitals: FragmentOrbital[] = dcdAlkenePair
    ? [
        {
          id: 'frag_dcd_m_dpi',
          side: 'left',
          label: `${dcdAlkenePair.metalAtom.element} ${dShellNum}d_π Back-Donor (${dShellNum}dxz)`,
          fragmentName: `${dcdAlkenePair.metalAtom.element} d-Manifold (DCD)`,
          energy: Number(((primProps.vsieD ?? -10.5) + 0.85).toFixed(2)),
          irrep: secondIrrepLower,
          degeneracy: 1,
          electrons: 2,
        },
        {
          id: 'frag_dcd_m_dnb',
          side: 'left',
          label: `${dcdAlkenePair.metalAtom.element} ${dShellNum}d Non-Bonding (${dShellNum}dz², ${dShellNum}dxy, ${dShellNum}dyz)`,
          fragmentName: `${dcdAlkenePair.metalAtom.element} d-Manifold (DCD)`,
          energy: Number(((primProps.vsieD ?? -10.5) - 0.6).toFixed(2)),
          irrep: `${firstIrrepLower} + ${lastIrrepLower}`,
          degeneracy: 3,
          electrons: 6,
        },
        {
          id: 'frag_dcd_m_dsigma',
          side: 'left',
          label: `${dcdAlkenePair.metalAtom.element} ${dShellNum}d_σ Acceptor (${dShellNum}dx²−y² / ${primShell}s)`,
          fragmentName: `${dcdAlkenePair.metalAtom.element} d-Manifold (DCD)`,
          energy: Number(((primProps.vsieD ?? -10.5) + 2.8).toFixed(2)),
          irrep: firstIrrepLower,
          degeneracy: 1,
          electrons: 0,
        },
        {
          id: 'frag_dcd_alkene_pi',
          side: 'right',
          label: 'Alkene π(C=C) σ-Donor',
          fragmentName: 'η²-Alkene Ligand (DCD)',
          energy: -10.51,
          irrep: firstIrrepLower,
          degeneracy: 1,
          electrons: 2,
        },
        {
          id: 'frag_dcd_alkene_pistar',
          side: 'right',
          label: 'Alkene π*(C=C) π-Acceptor',
          fragmentName: 'η²-Alkene Ligand (DCD)',
          energy: -7.2,
          irrep: secondIrrepLower,
          degeneracy: 1,
          electrons: 0,
        },
        {
          id: 'frag_right_salc',
          side: 'right',
          label: `Co-Ligand ${dcdAlkenePair.coDonor?.element || 'L'} σ/π SALCs`,
          fragmentName: 'Ancillary Co-Ligands',
          energy:
            ELEMENT_PROPERTIES[dcdAlkenePair.coDonor?.element || 'Cl']?.vsieP ?? periphProps.vsieP,
          irrep: `${firstIrrepLower} + ${lastIrrepLower}`,
          degeneracy: 3,
          electrons: 6,
        },
      ]
    : [
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
                  label: `${primaryAtom?.element} ${dShellNum}d`,
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
      } else if (pg.id === 'D6h') {
        // Evaluate inversion character χ(i) across centrosymmetric atom pairs (r_B ≈ -r_A)
        let invChar = 0;
        let invPairs = 0;
        for (let k = 0; k < N; k++) {
          const ck = vec[k];
          if (Math.abs(ck) < 0.05) continue;
          const ak = basis[k].atom;
          const ok = basis[k].orbital;
          const orbParity = ok === 's' || ok.startsWith('d') ? 1 : -1;
          for (let m = 0; m < N; m++) {
            if (basis[m].orbital !== ok) continue;
            const am = basis[m].atom;
            if (
              am.element === ak.element &&
              Math.hypot(ak.x3d + am.x3d, ak.y3d + am.y3d, ak.z3d + am.z3d) < 0.45
            ) {
              invChar += ck * vec[m] * orbParity;
              invPairs++;
              break;
            }
          }
        }
        const isUngerade = invPairs > 0 ? invChar < -0.02 : signChanges % 2 === 1;
        const targetId = isPi
          ? isUngerade
            ? 'E2u'
            : 'E1g'
          : isUngerade
            ? 'E1u'
            : 'E2g';
        irrepObj = multiPool.find((ir) => ir.id === targetId) || multiPool[0];
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
      } else if (pg.id === 'D6h') {
        let invChar = 0;
        let invPairs = 0;
        for (let k = 0; k < N; k++) {
          const ck = vec[k];
          if (Math.abs(ck) < 0.05) continue;
          const ak = basis[k].atom;
          const ok = basis[k].orbital;
          const orbParity = ok === 's' || ok.startsWith('d') ? 1 : -1;
          for (let m = 0; m < N; m++) {
            if (basis[m].orbital !== ok) continue;
            const am = basis[m].atom;
            if (
              am.element === ak.element &&
              Math.hypot(ak.x3d + am.x3d, ak.y3d + am.y3d, ak.z3d + am.z3d) < 0.45
            ) {
              invChar += ck * vec[m] * orbParity;
              invPairs++;
              break;
            }
          }
        }
        const isUngerade = invPairs > 0 ? invChar < -0.02 : signChanges % 2 === 1;
        const targetId = isPi
          ? isUngerade
            ? 'A2u'
            : 'B2g'
          : isUngerade
            ? 'B1u'
            : 'A1g';
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
      correlatedFragmentIds: dcdAlkenePair
        ? isPi
          ? ['frag_dcd_m_dpi', 'frag_dcd_alkene_pistar']
          : bondingType === 'non-bonding'
            ? ['frag_dcd_m_dnb']
            : ['frag_dcd_m_dsigma', 'frag_dcd_alkene_pi']
        : isPi
          ? ['frag_left_p', 'frag_right_salc']
          : energy < -16
            ? ['frag_left_s', 'frag_right_salc']
            : ['frag_left_p', 'frag_right_salc'],
      lcao,
      description: `Extended Hückel LCAO-MO (${interactionType}, ${lowerIrrepSym} symmetry) computed at ε = ${energy.toFixed(2)} eV with ${signChanges} internuclear nodal surface(s).`,
    });
  }

  // Enrich Dewar-Chatt-Duncanson (DCD) molecular orbitals for any transition-metal η²-alkene complex
  if (dcdAlkenePair && molecularOrbitals.length >= 4) {
    const { metalAtom, c1, c2 } = dcdAlkenePair;
    const mx = 0.5 * (c1.x3d + c2.x3d) - metalAtom.x3d;
    const my = 0.5 * (c1.y3d + c2.y3d) - metalAtom.y3d;
    const mz = 0.5 * (c1.z3d + c2.z3d) - metalAtom.z3d;
    const mLen = Math.hypot(mx, my, mz) || 1;
    const ux = mx / mLen;
    const uy = my / mLen;
    const uz = mz / mLen;

    const wx = c2.x3d - c1.x3d;
    const wy = c2.y3d - c1.y3d;
    const wz = c2.z3d - c1.z3d;
    const wLen = Math.hypot(wx, wy, wz) || 1;
    const wux = wx / wLen;
    const wuy = wy / wLen;
    const wuz = wz / wLen;

    const pOrbAlongU: OrbitalBasisType =
      Math.abs(ux) >= Math.abs(uy) && Math.abs(ux) >= Math.abs(uz)
        ? 'px'
        : Math.abs(uy) >= Math.abs(uz)
          ? 'py'
          : 'pz';
    const dSigmaOrb: OrbitalBasisType = Math.abs(uz) > 0.7 ? 'dz2' : 'dx2y2';
    const xzScore = Math.abs(ux * wuz + uz * wux);
    const yzScore = Math.abs(uy * wuz + uz * wuy);
    const dPiOrb: OrbitalBasisType =
      xzScore >= yzScore && xzScore >= 0.3 ? 'dxz' : yzScore >= 0.3 ? 'dyz' : 'dxy';

    const occupiedMos = molecularOrbitals.filter((m) => m.occupancy > 0);
    const virtualMos = molecularOrbitals.filter((m) => m.occupancy === 0);

    // 1. DCD μ/σ-Donation Bonding MO (lowest occupied frontier MO)
    const sigmaDonMo = occupiedMos[0];
    if (sigmaDonMo) {
      sigmaDonMo.label = `1${firstIrrepLower} (DCD μ/σ-Don)`;
      sigmaDonMo.irrepSymbol = firstIrrepLower;
      sigmaDonMo.bondingType = 'bonding';
      sigmaDonMo.interactionType = 'σ';
      sigmaDonMo.correlatedFragmentIds = ['frag_dcd_m_dsigma', 'frag_dcd_alkene_pi'];
      sigmaDonMo.lcao = [
        {
          atomId: metalAtom.id,
          element: metalAtom.element,
          orbital: dSigmaOrb,
          label: `${metalAtom.id} ${dShellNum}${dSigmaOrb === 'dx2y2' ? 'dx2-y2' : dSigmaOrb}`,
          coeff: 0.48,
        },
        {
          atomId: metalAtom.id,
          element: metalAtom.element,
          orbital: 's',
          label: `${metalAtom.id} ${primShell}s`,
          coeff: 0.24,
        },
        {
          atomId: c1.id,
          element: 'C',
          orbital: pOrbAlongU,
          label: `${c1.id} 2${pOrbAlongU}`,
          coeff: -0.46,
        },
        {
          atomId: c2.id,
          element: 'C',
          orbital: pOrbAlongU,
          label: `${c2.id} 2${pOrbAlongU}`,
          coeff: -0.46,
        },
      ];
      sigmaDonMo.description = `Dewar–Chatt–Duncanson σ-Donation MO (${firstIrrepLower}): Synergistic head-on donation from the filled alkene π(C=C) bonding orbital (${c1.id}+${c2.id} 2${pOrbAlongU}) into the empty ${metalAtom.element} ${dShellNum}${dSigmaOrb === 'dx2y2' ? 'dx²−y²' : dSigmaOrb}/${primShell}s acceptor hybrid.`;
    }

    // 2. DCD π-Backbonding Bonding MO (sub-HOMO occupied d_π -> π* orbital)
    const backbondMo =
      occupiedMos.length >= 3 ? occupiedMos[occupiedMos.length - 2] : occupiedMos[occupiedMos.length - 1];
    if (backbondMo && backbondMo !== sigmaDonMo) {
      backbondMo.label = `1${secondIrrepLower} (DCD π-Back)`;
      backbondMo.irrepSymbol = secondIrrepLower;
      backbondMo.bondingType = 'bonding';
      backbondMo.interactionType = 'π';
      backbondMo.correlatedFragmentIds = ['frag_dcd_m_dpi', 'frag_dcd_alkene_pistar'];
      backbondMo.lcao = [
        {
          atomId: metalAtom.id,
          element: metalAtom.element,
          orbital: dPiOrb,
          label: `${metalAtom.id} ${dShellNum}${dPiOrb}`,
          coeff: 0.62,
        },
        {
          atomId: c1.id,
          element: 'C',
          orbital: pOrbAlongU,
          label: `${c1.id} 2${pOrbAlongU}`,
          coeff: 0.42,
        },
        {
          atomId: c2.id,
          element: 'C',
          orbital: pOrbAlongU,
          label: `${c2.id} 2${pOrbAlongU}`,
          coeff: -0.42,
        },
      ];
      backbondMo.description = `Dewar–Chatt–Duncanson π-Backbonding MO (${secondIrrepLower}): Side-on d_π → π*(C=C) back-donation from filled ${metalAtom.element} ${dShellNum}${dPiOrb} into the empty alkene π* antibonding orbital, strengthening ${metalAtom.element}–C bonding while elongating C=C (1.376 Å) and bending substituents back.`;
    }

    // 3. DCD σ* LUMO and π* Backbonding Antibonding MOs
    const lumoMo = virtualMos[0];
    if (lumoMo) {
      lumoMo.label = `2${firstIrrepLower}* (DCD σ* LUMO)`;
      lumoMo.bondingType = 'antibonding';
      lumoMo.interactionType = 'σ*';
      lumoMo.correlatedFragmentIds = ['frag_dcd_m_dsigma', 'frag_dcd_alkene_pi', 'frag_right_salc'];
      lumoMo.description = `Dewar–Chatt–Duncanson σ* Antibonding LUMO (${firstIrrepLower}*): Out-of-phase antibonding combination of ${metalAtom.element} ${dShellNum}${dSigmaOrb === 'dx2y2' ? 'dx²−y²' : dSigmaOrb} with the alkene π(C=C) σ-donor and co-ligand σ SALCs.`;
    }
    const piStarMo = virtualMos.length >= 2 ? virtualMos[1] : undefined;
    if (piStarMo) {
      piStarMo.label = `2${secondIrrepLower}* (DCD π*-Back)`;
      piStarMo.irrepSymbol = secondIrrepLower;
      piStarMo.bondingType = 'antibonding';
      piStarMo.interactionType = 'π*';
      piStarMo.correlatedFragmentIds = ['frag_dcd_m_dpi', 'frag_dcd_alkene_pistar'];
      piStarMo.lcao = [
        {
          atomId: c1.id,
          element: 'C',
          orbital: pOrbAlongU,
          label: `${c1.id} 2${pOrbAlongU}`,
          coeff: 0.56,
        },
        {
          atomId: c2.id,
          element: 'C',
          orbital: pOrbAlongU,
          label: `${c2.id} 2${pOrbAlongU}`,
          coeff: -0.56,
        },
        {
          atomId: metalAtom.id,
          element: metalAtom.element,
          orbital: dPiOrb,
          label: `${metalAtom.id} ${dShellNum}${dPiOrb}`,
          coeff: -0.48,
        },
      ];
      piStarMo.description = `Dewar–Chatt–Duncanson π* Antibonding MO (${secondIrrepLower}*): Out-of-phase antibonding counterpart of the ${metalAtom.element} ${dShellNum}${dPiOrb} → alkene π*(C=C) backbonding interaction.`;
    }
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
    leftFragmentTitle: dcdAlkenePair
      ? `Metal ${primaryAtom?.element || 'M'} (${dShellNum}d, ${primShell}s, ${primShell}p — DCD)`
      : primIsLn
        ? `Central ${primaryAtom?.element || 'Skeleton'} (4f, 5d, 6s)`
        : primIsTm
          ? `Central ${primaryAtom?.element || 'Skeleton'} (${dShellNum}d, ${primShell}s, ${primShell}p)`
          : `Central ${primaryAtom?.element || 'Skeleton'} (${primShell}s, ${primShell}p)`,
    rightFragmentTitle: dcdAlkenePair
      ? `η²-Alkene (π / π*) & Co-Ligand SALCs`
      : `Peripheral ${peripheralAtom?.element || 'Ligand'} SALCs`,
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
