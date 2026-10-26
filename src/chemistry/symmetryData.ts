export interface SymmetryClass {
  id: string;
  symbol: string; // e.g., "E", "2C3", "3σv"
  operationType: 'E' | 'Cn' | 'Sn' | 'sigma_h' | 'sigma_v' | 'sigma_d' | 'i';
  count: number; // N_k in character table
  description: string;
}

export interface IrreducibleRepresentation {
  symbol: string; // Mulliken symbol e.g. "A1", "B1", "E", "T2", "E1g"
  id: string;     // ASCII key e.g. "A1", "B1", "E", "T2", "E1g"
  dimension: number; // 1 for A/B, 2 for E, 3 for T
  characters: number[]; // Matches pointGroup.classes length
  linearFunctions: string[]; // e.g. ["z"], ["x", "Ry"]
  quadraticFunctions: string[]; // e.g. ["x2", "y2", "z2"], ["xy"]
  orbitalBasis: string; // e.g. "s, pz, dz2"
  irActive: boolean;
  ramanActive: boolean;
  description: string;
}

export interface SymmetryElement3D {
  id: string;
  type: 'axis' | 'plane' | 'center';
  label: string;
  schoenflies: string;
  fold?: number; // e.g. 2, 3, 4, 6 for Cn
  axis?: [number, number, number]; // unit vector for axis or plane normal
  color: string;
  description: string;
}

export interface ReduciblePreset {
  id: string;
  label: string;
  characters: number[];
  description: string;
}

export interface PointGroupData {
  id: string;
  schoenflies: string; // e.g. "C2v"
  name: string;
  order: number; // h
  family: string;
  description: string;
  classes: SymmetryClass[];
  irreps: IrreducibleRepresentation[];
  elements3D: SymmetryElement3D[];
  defaultReduciblePresets: ReduciblePreset[];
}

export interface DecompositionResult {
  irrepSymbol: string;
  irrepId: string;
  coefficient: number;
  sumContribution: number;
  terms: {
    classSymbol: string;
    count: number;
    reducibleChar: number;
    irrepChar: number;
    product: number;
  }[];
}

export const POINT_GROUPS: Record<string, PointGroupData> = {
  C2v: {
    id: 'C2v',
    schoenflies: 'C2v',
    name: 'Dihedral Pyramidal / Planar Bent',
    order: 4,
    family: 'Cnv (Pyramidal / Polar)',
    description: 'Contains a principal C2 rotation axis and two mutually perpendicular vertical mirror planes σv(xz) and σv\'(yz). Classic examples: H2O, H2CO, SO2, CH2Cl2.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation (leaves all coordinates unchanged)' },
      { id: 'C2', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° proper rotation around the principal z-axis' },
      { id: 'sigma_xz', symbol: 'σv(xz)', operationType: 'sigma_v', count: 1, description: 'Vertical reflection across the xz molecular plane' },
      { id: 'sigma_yz', symbol: 'σv\'(yz)', operationType: 'sigma_v', count: 1, description: 'Vertical reflection across the yz plane bisecting the bond angle' },
    ],
    irreps: [
      {
        symbol: 'A1',
        id: 'A1',
        dimension: 1,
        characters: [1, 1, 1, 1],
        linearFunctions: ['z'],
        quadraticFunctions: ['x²', 'y²', 'z²'],
        orbitalBasis: 's, pz, dz², dx²-y²',
        irActive: true,
        ramanActive: true,
        description: 'Totally symmetric representation; invariant under C2 rotation and both mirror planes.',
      },
      {
        symbol: 'A2',
        id: 'A2',
        dimension: 1,
        characters: [1, 1, -1, -1],
        linearFunctions: ['Rz'],
        quadraticFunctions: ['xy'],
        orbitalBasis: 'dxy',
        irActive: false,
        ramanActive: true,
        description: 'Symmetric with respect to C2 rotation, antisymmetric with respect to both vertical mirror planes.',
      },
      {
        symbol: 'B1',
        id: 'B1',
        dimension: 1,
        characters: [1, -1, 1, -1],
        linearFunctions: ['x', 'Ry'],
        quadraticFunctions: ['xz'],
        orbitalBasis: 'px, dxz',
        irActive: true,
        ramanActive: true,
        description: 'Antisymmetric under C2 rotation, symmetric across σv(xz), antisymmetric across σv\'(yz).',
      },
      {
        symbol: 'B2',
        id: 'B2',
        dimension: 1,
        characters: [1, -1, -1, 1],
        linearFunctions: ['y', 'Rx'],
        quadraticFunctions: ['yz'],
        orbitalBasis: 'py, dyz',
        irActive: true,
        ramanActive: true,
        description: 'Antisymmetric under C2 rotation and σv(xz), symmetric in the molecular yz plane.',
      },
    ],
    elements3D: [
      { id: 'c2_z', type: 'axis', label: 'C2 (z-axis)', schoenflies: 'C2', fold: 2, axis: [0, 1, 0], color: '#2563EB', description: '180° rotational symmetry axis along the molecular bisector' },
      { id: 'sigma_xz', type: 'plane', label: 'σv(xz)', schoenflies: 'σv', axis: [0, 0, 1], color: '#0D9488', description: 'Mirror plane perpendicular to the molecular plane' },
      { id: 'sigma_yz', type: 'plane', label: 'σv\'(yz)', schoenflies: 'σv\'', axis: [1, 0, 0], color: '#D97706', description: 'Mirror plane containing all atoms in planar C2v molecules' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_bonds', label: 'Γ_σ (2 O–H σ bonds)', characters: [2, 0, 0, 2], description: 'Basis of the two σ bond vectors in a bent AB2 molecule (Γ_σ = A1 + B2)' },
      { id: 'valence_ao', label: 'Γ_val (H2O 6 Valence AOs)', characters: [6, 0, 2, 4], description: 'Oxygen (2s, 2px, 2py, 2pz) + two Hydrogen (1s) atomic orbitals (Γ_val = 3A1 + B1 + 2B2)' },
      { id: 'cartesian_3N', label: 'Γ_3N (9 Cartesian Coordinates)', characters: [9, -1, 1, 3], description: 'Full 3N degrees of freedom for triatomic AB2 (Γ_3N = 3A1 + A2 + 2B1 + 3B2)' },
    ],
  },

  C3v: {
    id: 'C3v',
    schoenflies: 'C3v',
    name: 'Trigonal Pyramidal',
    order: 6,
    family: 'Cnv (Pyramidal / Polar)',
    description: 'Features a 3-fold principal axis C3 and three vertical mirror planes σv intersecting along the C3 axis. Classic examples: NH3, PH3, CHCl3, POCl3.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C3', symbol: '2C3(z)', operationType: 'Cn', count: 2, description: '±120° rotations around the principal z-axis' },
      { id: '3sigma_v', symbol: '3σv', operationType: 'sigma_v', count: 3, description: 'Three vertical mirror planes each containing the central atom and one ligand' },
    ],
    irreps: [
      {
        symbol: 'A1',
        id: 'A1',
        dimension: 1,
        characters: [1, 1, 1],
        linearFunctions: ['z'],
        quadraticFunctions: ['x²+y²', 'z²'],
        orbitalBasis: 's, pz, dz²',
        irActive: true,
        ramanActive: true,
        description: 'Totally symmetric representation along the C3 axis.',
      },
      {
        symbol: 'A2',
        id: 'A2',
        dimension: 1,
        characters: [1, 1, -1],
        linearFunctions: ['Rz'],
        quadraticFunctions: [],
        orbitalBasis: 'fx(x²-3y²)',
        irActive: false,
        ramanActive: false,
        description: 'Symmetric under C3 rotation, antisymmetric under vertical reflections σv.',
      },
      {
        symbol: 'E',
        id: 'E',
        dimension: 2,
        characters: [2, -1, 0],
        linearFunctions: ['(x, y)', '(Rx, Ry)'],
        quadraticFunctions: ['(x²-y², xy)', '(xz, yz)'],
        orbitalBasis: '(px, py), (dxz, dyz), (dx²-y², dxy)',
        irActive: true,
        ramanActive: true,
        description: 'Doubly degenerate representation spanning the xy plane perpendicular to C3.',
      },
    ],
    elements3D: [
      { id: 'c3_z', type: 'axis', label: 'C3 (z-axis)', schoenflies: 'C3', fold: 3, axis: [0, 1, 0], color: '#2563EB', description: '120° 3-fold principal rotation axis' },
      { id: 'sigma_v1', type: 'plane', label: 'σv(1)', schoenflies: 'σv', axis: [0, 0, 1], color: '#0D9488', description: 'First vertical reflection plane containing N–H1' },
      { id: 'sigma_v2', type: 'plane', label: 'σv(2)', schoenflies: 'σv', axis: [0.866, 0, -0.5], color: '#0D9488', description: 'Second vertical reflection plane rotated by 120°' },
      { id: 'sigma_v3', type: 'plane', label: 'σv(3)', schoenflies: 'σv', axis: [-0.866, 0, -0.5], color: '#0D9488', description: 'Third vertical reflection plane rotated by 240°' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_bonds', label: 'Γ_σ (3 N–H σ bonds)', characters: [3, 0, 1], description: 'Three N–H sigma bonds in NH3 decompose into Γ_σ = A1 + E' },
      { id: 'valence_ao', label: 'Γ_val (NH3 7 Valence AOs)', characters: [7, 1, 3], description: 'N(2s, 2p) + 3 H(1s) orbitals decompose into Γ_val = 3A1 + 2E' },
      { id: 'cartesian_3N', label: 'Γ_3N (12 Cartesian Coordinates)', characters: [12, 0, 2], description: 'Full 3N coordinates for pyramidal AB3 (Γ_3N = 3A1 + A2 + 4E)' },
    ],
  },

  D2h: {
    id: 'D2h',
    schoenflies: 'D2h',
    name: 'Orthorhombic Prismatic / Planar Rectangular',
    order: 8,
    family: 'Dnh (Centrosymmetric Dihedral)',
    description: 'Three mutually perpendicular C2 axes, an inversion center i, and three mutually perpendicular mirror planes (xy, xz, yz). Classic examples: Ethylene (C2H4), Diborane (B2H6), N2O4.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: 'C2z', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation about the C=C bond axis (z)' },
      { id: 'C2y', symbol: 'C2(y)', operationType: 'Cn', count: 1, description: '180° rotation about the in-plane perpendicular y-axis' },
      { id: 'C2x', symbol: 'C2(x)', operationType: 'Cn', count: 1, description: '180° rotation about the out-of-plane x-axis' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the molecular centroid (x,y,z) → (-x,-y,-z)' },
      { id: 'sigma_xy', symbol: 'σ(xy)', operationType: 'sigma_h', count: 1, description: 'Reflection across the xy plane bisecting the C=C bond' },
      { id: 'sigma_xz', symbol: 'σ(xz)', operationType: 'sigma_v', count: 1, description: 'Reflection across the xz plane perpendicular to the molecular plane' },
      { id: 'sigma_yz', symbol: 'σ(yz)', operationType: 'sigma_v', count: 1, description: 'Reflection in the yz molecular plane' },
    ],
    irreps: [
      { symbol: 'Ag', id: 'Ag', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²', 'y²', 'z²'], orbitalBasis: 's, dz², dx²-y²', irActive: false, ramanActive: true, description: 'Totally symmetric gerade representation.' },
      { symbol: 'B1g', id: 'B1g', dimension: 1, characters: [1, 1, -1, -1, 1, 1, -1, -1], linearFunctions: ['Rz'], quadraticFunctions: ['xy'], orbitalBasis: 'dxy', irActive: false, ramanActive: true, description: 'Gerade representation symmetric under C2(z).' },
      { symbol: 'B2g', id: 'B2g', dimension: 1, characters: [1, -1, 1, -1, 1, -1, 1, -1], linearFunctions: ['Ry'], quadraticFunctions: ['xz'], orbitalBasis: 'dxz (π* antibonding)', irActive: false, ramanActive: true, description: 'Gerade representation corresponding to out-of-plane π* orbital in ethylene.' },
      { symbol: 'B3g', id: 'B3g', dimension: 1, characters: [1, -1, -1, 1, 1, -1, -1, 1], linearFunctions: ['Rx'], quadraticFunctions: ['yz'], orbitalBasis: 'dyz', irActive: false, ramanActive: true, description: 'Gerade representation symmetric in the yz molecular plane.' },
      { symbol: 'Au', id: 'Au', dimension: 1, characters: [1, 1, 1, 1, -1, -1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'fxyz', irActive: false, ramanActive: false, description: 'Ungerade twisting representation (silent in both IR and Raman).' },
      { symbol: 'B1u', id: 'B1u', dimension: 1, characters: [1, 1, -1, -1, -1, -1, 1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz', irActive: true, ramanActive: false, description: 'Ungerade representation along the principal C–C axis.' },
      { symbol: 'B2u', id: 'B2u', dimension: 1, characters: [1, -1, 1, -1, -1, 1, -1, 1], linearFunctions: ['y'], quadraticFunctions: [], orbitalBasis: 'py', irActive: true, ramanActive: false, description: 'Ungerade in-plane representation.' },
      { symbol: 'B3u', id: 'B3u', dimension: 1, characters: [1, -1, -1, 1, -1, 1, 1, -1], linearFunctions: ['x'], quadraticFunctions: [], orbitalBasis: 'px (π bonding)', irActive: true, ramanActive: false, description: 'Ungerade out-of-plane π bonding representation in ethylene.' },
    ],
    elements3D: [
      { id: 'c2_z', type: 'axis', label: 'C2(z)', schoenflies: 'C2(z)', fold: 2, axis: [1, 0, 0], color: '#2563EB', description: 'C2 axis along the C=C bond' },
      { id: 'c2_y', type: 'axis', label: 'C2(y)', schoenflies: 'C2(y)', fold: 2, axis: [0, 1, 0], color: '#2563EB', description: 'In-plane perpendicular C2 axis' },
      { id: 'c2_x', type: 'axis', label: 'C2(x)', schoenflies: 'C2(x)', fold: 2, axis: [0, 0, 1], color: '#2563EB', description: 'Out-of-plane perpendicular C2 axis' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion at midpoint of C=C bond' },
      { id: 'sigma_yz', type: 'plane', label: 'σ(yz)', schoenflies: 'σ(yz)', axis: [0, 0, 1], color: '#0D9488', description: 'Molecular plane containing all 6 atoms' },
      { id: 'sigma_xy', type: 'plane', label: 'σ(xy)', schoenflies: 'σ(xy)', axis: [1, 0, 0], color: '#D97706', description: 'Mirror plane bisecting the C=C bond' },
    ],
    defaultReduciblePresets: [
      { id: 'pi_system', label: 'Γ_π (2 Carbon 2px π Orbitals)', characters: [2, -2, 0, 0, 0, 2, 0, -2], description: 'Two out-of-plane 2px orbitals form bonding b3u (π) and antibonding b2g (π*)' },
      { id: 'ch_bonds', label: 'Γ_CH (4 C–H σ Bonds)', characters: [4, 0, 0, 0, 0, 0, 0, 4], description: 'Four C–H sigma bonds decompose into Ag + B3g + B1u + B2u' },
    ],
  },

  D3h: {
    id: 'D3h',
    schoenflies: 'D3h',
    name: 'Trigonal Planar / Trigonal Bipyramidal',
    order: 12,
    family: 'Dnh (Dihedral with Horizontal Mirror)',
    description: 'Contains a 3-fold principal C3 axis, three perpendicular C2 axes, a horizontal mirror plane σh, a 3-fold improper axis S3, and three vertical mirror planes σv. Examples: BF3, BH3, CO3(2-), PCl5.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C3', symbol: '2C3(z)', operationType: 'Cn', count: 2, description: '±120° rotations around the out-of-plane z-axis' },
      { id: '3C2', symbol: '3C2\'', operationType: 'Cn', count: 3, description: 'Three 180° rotations along each B–F / B–H bond axis' },
      { id: 'sigma_h', symbol: 'σh(xy)', operationType: 'sigma_h', count: 1, description: 'Horizontal reflection in the trigonal molecular plane' },
      { id: '2S3', symbol: '2S3', operationType: 'Sn', count: 2, description: '120° rotation followed by horizontal reflection σh' },
      { id: '3sigma_v', symbol: '3σv', operationType: 'sigma_v', count: 3, description: 'Three vertical mirror planes containing each bond axis' },
    ],
    irreps: [
      { symbol: 'A1\'', id: 'A1p', dimension: 1, characters: [1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation in the trigonal plane.' },
      { symbol: 'A2\'', id: 'A2p', dimension: 1, characters: [1, 1, -1, 1, 1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'In-plane rotational symmetry.' },
      { symbol: 'E\'', id: 'Ep', dimension: 2, characters: [2, -1, 0, 2, -1, 0], linearFunctions: ['(x, y)'], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(px, py), (dx²-y², dxy)', irActive: true, ramanActive: true, description: 'Doubly degenerate in-plane σ bonding representation.' },
      { symbol: 'A1\'\'', id: 'A1pp', dimension: 1, characters: [1, 1, 1, -1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Antisymmetric across σh and σv.' },
      { symbol: 'A2\'\'', id: 'A2pp', dimension: 1, characters: [1, 1, -1, -1, -1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz (out-of-plane π)', irActive: true, ramanActive: false, description: 'Out-of-plane pz orbital (empty LUMO on BH3 or π bonding in BF3).' },
      { symbol: 'E\'\'', id: 'Epp', dimension: 2, characters: [2, -1, 0, -2, 1, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz)', irActive: false, ramanActive: true, description: 'Doubly degenerate out-of-plane π representation.' },
    ],
    elements3D: [
      { id: 'c3_z', type: 'axis', label: 'C3 (z-axis)', schoenflies: 'C3', fold: 3, axis: [0, 0, 1], color: '#2563EB', description: 'Out-of-plane 3-fold rotation axis' },
      { id: 'c2_1', type: 'axis', label: 'C2\'(1)', schoenflies: 'C2\'', fold: 2, axis: [0, 1, 0], color: '#0284C7', description: 'In-plane C2 axis along bond 1' },
      { id: 'sigma_h', type: 'plane', label: 'σh (xy plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#D97706', description: 'Horizontal mirror plane containing all atoms' },
      { id: 'sigma_v1', type: 'plane', label: 'σv(1)', schoenflies: 'σv', axis: [1, 0, 0], color: '#0D9488', description: 'Vertical mirror plane along bond 1' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_bonds', label: 'Γ_σ (3 Trigonal σ Bonds)', characters: [3, 0, 1, 3, 0, 1], description: 'Three in-plane σ bonds decompose into Γ_σ = A1\' + E\' (sp² hybridization)' },
      { id: 'pi_bonds', label: 'Γ_π (3 Out-of-plane pz Orbitals)', characters: [3, 0, -1, -3, 0, 1], description: 'Three ligand pz orbitals decompose into Γ_π = A2\'\' + E\'\'' },
    ],
  },

  D4h: {
    id: 'D4h',
    schoenflies: 'D4h',
    name: 'Square Planar / Octahedral Trans-Substituted',
    order: 16,
    family: 'Dnh (Centrosymmetric Dihedral)',
    description: '4-fold principal C4 axis, four perpendicular C2 axes (2 C2\' along bonds, 2 C2\'\' bisecting bonds), inversion center i, and horizontal mirror plane σh. Classic examples: XeF4, [PtCl4]2-, Cyclobutadiene.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C4', symbol: '2C4(z)', operationType: 'Cn', count: 2, description: '±90° rotations around the z-axis' },
      { id: 'C2', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation around the z-axis' },
      { id: '2C2p', symbol: '2C2\'', operationType: 'Cn', count: 2, description: 'Two 180° rotations along the ligand axes (x and y)' },
      { id: '2C2pp', symbol: '2C2\'\'', operationType: 'Cn', count: 2, description: 'Two 180° rotations bisecting the ligand axes' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Center of inversion' },
      { id: '2S4', symbol: '2S4', operationType: 'Sn', count: 2, description: '90° improper rotations' },
      { id: 'sigma_h', symbol: 'σh', operationType: 'sigma_h', count: 1, description: 'Horizontal mirror plane (xy)' },
      { id: '2sigma_v', symbol: '2σv', operationType: 'sigma_v', count: 2, description: 'Vertical mirror planes containing bond axes' },
      { id: '2sigma_d', symbol: '2σd', operationType: 'sigma_d', count: 2, description: 'Dihedral mirror planes bisecting bond axes' },
    ],
    irreps: [
      { symbol: 'A1g', id: 'A1g', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation.' },
      { symbol: 'B1g', id: 'B1g', dimension: 1, characters: [1, -1, 1, 1, -1, 1, -1, 1, 1, -1], linearFunctions: [], quadraticFunctions: ['x²-y²'], orbitalBasis: 'dx²-y²', irActive: false, ramanActive: true, description: 'Matches lobes pointing directly along the x and y ligand axes.' },
      { symbol: 'B2g', id: 'B2g', dimension: 1, characters: [1, -1, 1, -1, 1, 1, -1, 1, -1, 1], linearFunctions: [], quadraticFunctions: ['xy'], orbitalBasis: 'dxy', irActive: false, ramanActive: true, description: 'In-plane dxy orbital bisecting the bonds.' },
      { symbol: 'Eg', id: 'Eg', dimension: 2, characters: [2, 0, -2, 0, 0, 2, 0, -2, 0, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz)', irActive: false, ramanActive: true, description: 'Doubly degenerate out-of-plane d-orbital pair.' },
      { symbol: 'A2u', id: 'A2u', dimension: 1, characters: [1, 1, 1, -1, -1, -1, -1, -1, 1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz', irActive: true, ramanActive: false, description: 'Out-of-plane ungerade pz orbital.' },
      { symbol: 'Eu', id: 'Eu', dimension: 2, characters: [2, 0, -2, 0, 0, -2, 0, 2, 0, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py)', irActive: true, ramanActive: false, description: 'Doubly degenerate in-plane px, py orbitals.' },
    ],
    elements3D: [
      { id: 'c4_z', type: 'axis', label: 'C4 (z-axis)', schoenflies: 'C4', fold: 4, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 4-fold axis perpendicular to square plane' },
      { id: 'sigma_h', type: 'plane', label: 'σh (xy plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#D97706', description: 'Square planar mirror plane' },
      { id: 'sigma_v1', type: 'plane', label: 'σv(xz)', schoenflies: 'σv', axis: [0, 1, 0], color: '#0D9488', description: 'Vertical mirror plane along x-axis ligands' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Central atom inversion center' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_bonds', label: 'Γ_σ (4 Square Planar σ Bonds)', characters: [4, 0, 0, 2, 0, 0, 0, 4, 2, 0], description: 'Four σ bonds decompose into Γ_σ = A1g + B1g + Eu (dsp² hybridization)' },
    ],
  },

  D6h: {
    id: 'D6h',
    schoenflies: 'D6h',
    name: 'Hexagonal Planar (Aromatic)',
    order: 24,
    family: 'Dnh (Centrosymmetric Dihedral)',
    description: '6-fold principal C6 axis perpendicular to the ring, six perpendicular C2 axes (3 C2\' through opposite atoms, 3 C2\'\' bisecting opposite bonds), inversion center i, improper axes 2S3 and 2S6, and horizontal plane σh. Classic example: Benzene (C6H6).',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C6', symbol: '2C6(z)', operationType: 'Cn', count: 2, description: '±60° rotations around the ring normal' },
      { id: '2C3', symbol: '2C3(z)', operationType: 'Cn', count: 2, description: '±120° rotations around the ring normal' },
      { id: 'C2', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation around the ring normal' },
      { id: '3C2p', symbol: '3C2\'', operationType: 'Cn', count: 3, description: 'In-plane 180° axes passing through opposite C–H vertices' },
      { id: '3C2pp', symbol: '3C2\'\'', operationType: 'Cn', count: 3, description: 'In-plane 180° axes bisecting opposite C–C bonds' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the center of the ring' },
      { id: '2S3', symbol: '2S3', operationType: 'Sn', count: 2, description: '±120° improper rotations around z-axis' },
      { id: '2S6', symbol: '2S6', operationType: 'Sn', count: 2, description: '±60° improper rotations around z-axis' },
      { id: 'sigma_h', symbol: 'σh', operationType: 'sigma_h', count: 1, description: 'Horizontal mirror plane containing the hexagonal ring' },
      { id: '3sigma_d', symbol: '3σd', operationType: 'sigma_d', count: 3, description: 'Dihedral mirror planes bisecting opposite C–C bonds' },
      { id: '3sigma_v', symbol: '3σv', operationType: 'sigma_v', count: 3, description: 'Vertical mirror planes through opposite C–H vertices' },
    ],
    irreps: [
      { symbol: 'A1g', id: 'A1g', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric ring breathing & lowest σ orbital.' },
      { symbol: 'B2g', id: 'B2g', dimension: 1, characters: [1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'fz(x²-3y²) (Highest π* MO)', irActive: false, ramanActive: false, description: 'Highest energy 3-node π* antibonding orbital of benzene (b2g).' },
      { symbol: 'E1g', id: 'E1g', dimension: 2, characters: [2, 1, -1, -2, 0, 0, 2, 1, -1, -2, 0, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz) — Degenerate π HOMO', irActive: false, ramanActive: true, description: 'Doubly degenerate 1-node π bonding HOMO of benzene (e1g).' },
      { symbol: 'E2g', id: 'E2g', dimension: 2, characters: [2, -1, -1, 2, 0, 0, 2, -1, -1, 2, 0, 0], linearFunctions: [], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(dx²-y², dxy)', irActive: false, ramanActive: true, description: 'Doubly degenerate in-plane ring deformation.' },
      { symbol: 'A2u', id: 'A2u', dimension: 1, characters: [1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz — Lowest π MO', irActive: true, ramanActive: false, description: 'Nodeless out-of-plane π bonding MO of benzene (a2u).' },
      { symbol: 'B1u', id: 'B1u', dimension: 1, characters: [1, -1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'In-plane C–C σ', irActive: false, ramanActive: false, description: 'Antisymmetric ungerade in-plane ring mode.' },
      { symbol: 'E1u', id: 'E1u', dimension: 2, characters: [2, 1, -1, -2, 0, 0, -2, -1, 1, 2, 0, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py)', irActive: true, ramanActive: false, description: 'Doubly degenerate in-plane IR-active ring modes.' },
      { symbol: 'E2u', id: 'E2u', dimension: 2, characters: [2, -1, -1, 2, 0, 0, -2, 1, 1, -2, 0, 0], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'Degenerate π* LUMO', irActive: false, ramanActive: false, description: 'Doubly degenerate 2-node π* antibonding LUMO of benzene (e2u).' },
    ],
    elements3D: [
      { id: 'c6_z', type: 'axis', label: 'C6 (z-axis)', schoenflies: 'C6', fold: 6, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 6-fold axis perpendicular to the benzene ring' },
      { id: 'c2_p1', type: 'axis', label: 'C2\' (vertex axis)', schoenflies: 'C2\'', fold: 2, axis: [1, 0, 0], color: '#0284C7', description: 'In-plane 2-fold axis through opposite carbons' },
      { id: 'sigma_h', type: 'plane', label: 'σh (Ring Plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#D97706', description: 'Horizontal mirror plane containing all 12 atoms (nodal plane of all π MOs)' },
      { id: 'sigma_v1', type: 'plane', label: 'σv (Vertex Plane)', schoenflies: 'σv', axis: [0, 1, 0], color: '#0D9488', description: 'Vertical mirror plane bisecting opposite vertices' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Ring centroid inversion center' },
    ],
    defaultReduciblePresets: [
      { id: 'pi_6', label: 'Γ_π (6 Carbon 2pz π Orbitals)', characters: [6, 0, 0, 0, -2, 0, 0, 0, 0, -6, 0, 2], description: 'Six 2pz orbitals decompose into Γ_π = B2g + E1g + A2u + E2u' },
      { id: 'sigma_cc', label: 'Γ_CC (6 C–C Ring σ Bonds)', characters: [6, 0, 0, 0, 0, 2, 0, 0, 0, 6, 2, 0], description: 'Six C–C ring σ bonds decompose into A1g + E2g + B1u + E1u' },
    ],
  },

  Td: {
    id: 'Td',
    schoenflies: 'Td',
    name: 'Tetrahedral',
    order: 24,
    family: 'Cubic / Polyhedral',
    description: 'Full symmetry of a regular tetrahedron: four C3 axes along each bond, three C2 axes bisecting bond angles, six S4 improper axes, and six dihedral mirror planes σd. Classic examples: CH4, CCl4, NH4+, SiH4.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '8C3', symbol: '8C3', operationType: 'Cn', count: 8, description: 'Eight 120° rotations along the four C–H / M–L bond axes' },
      { id: '3C2', symbol: '3C2', operationType: 'Cn', count: 3, description: 'Three 180° rotations along the Cartesian x, y, z axes' },
      { id: '6S4', symbol: '6S4', operationType: 'Sn', count: 6, description: 'Six 90° rotary-reflections along the x, y, z axes' },
      { id: '6sigma_d', symbol: '6σd', operationType: 'sigma_d', count: 6, description: 'Six dihedral mirror planes each containing two bonds' },
    ],
    irreps: [
      { symbol: 'A1', id: 'A1', dimension: 1, characters: [1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²+z²'], orbitalBasis: 's', irActive: false, ramanActive: true, description: 'Totally symmetric spherical representation (2s orbital on C).' },
      { symbol: 'A2', id: 'A2', dimension: 1, characters: [1, 1, 1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'fxyz', irActive: false, ramanActive: false, description: 'Antisymmetric under S4 and σd.' },
      { symbol: 'E', id: 'E', dimension: 2, characters: [2, -1, 2, 0, 0], linearFunctions: [], quadraticFunctions: ['(2z²-x²-y², x²-y²)'], orbitalBasis: '(dz², dx²-y²)', irActive: false, ramanActive: true, description: 'Doubly degenerate d-orbital pair pointing between tetrahedral ligands.' },
      { symbol: 'T1', id: 'T1', dimension: 3, characters: [3, 0, -1, 1, -1], linearFunctions: ['(Rx, Ry, Rz)'], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Triply degenerate rotational representation.' },
      { symbol: 'T2', id: 'T2', dimension: 3, characters: [3, 0, -1, -1, 1], linearFunctions: ['(x, y, z)'], quadraticFunctions: ['(xy, xz, yz)'], orbitalBasis: '(px, py, pz), (dxy, dxz, dyz)', irActive: true, ramanActive: true, description: 'Triply degenerate p and d orbitals (sp³ / sd³ tetrahedral bonding).' },
    ],
    elements3D: [
      { id: 'c3_1', type: 'axis', label: 'C3 (Bond Axis)', schoenflies: 'C3', fold: 3, axis: [0.577, 0.577, 0.577], color: '#2563EB', description: '3-fold axis along a vertex of the tetrahedron' },
      { id: 'c2_z', type: 'axis', label: 'C2 / S4 (z-axis)', schoenflies: 'S4/C2', fold: 4, axis: [0, 1, 0], color: '#0284C7', description: 'Collinear C2 and S4 axis bisecting opposite H–C–H angles' },
      { id: 'sigma_d1', type: 'plane', label: 'σd (Dihedral Plane)', schoenflies: 'σd', axis: [0.707, 0, -0.707], color: '#0D9488', description: 'Dihedral mirror plane containing two C–H bonds and bisecting the other two' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_4', label: 'Γ_σ (4 Tetrahedral C–H Bonds)', characters: [4, 1, 0, 0, 2], description: 'Four tetrahedral σ bonds decompose into Γ_σ = A1 + T2 (confirming sp³ hybridization)' },
      { id: 'cartesian_3N', label: 'Γ_3N (15 Cartesian Coordinates)', characters: [15, 0, -1, -1, 3], description: 'Full 15 degrees of freedom for CH4 decompose into A1 + E + T1 + 3T2' },
    ],
  },

  Oh: {
    id: 'Oh',
    schoenflies: 'Oh',
    name: 'Octahedral',
    order: 48,
    family: 'Cubic / Polyhedral',
    description: 'Full octahedral symmetry with inversion center i, three 4-fold C4 axes along the Cartesian x, y, z ligand axes, four 3-fold C3 axes through octahedral faces, and nine mirror planes (3 σh + 6 σd). Classic examples: SF6, [Fe(CN)6]4-, Mo(CO)6.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '8C3', symbol: '8C3', operationType: 'Cn', count: 8, description: '120° rotations through opposite triangular faces' },
      { id: '6C2', symbol: '6C2\'', operationType: 'Cn', count: 6, description: '180° rotations bisecting opposite edges' },
      { id: '6C4', symbol: '6C4', operationType: 'Cn', count: 6, description: '90° rotations along the x, y, z ligand axes' },
      { id: '3C2z', symbol: '3C2(=C4²)', operationType: 'Cn', count: 3, description: '180° rotations along the x, y, z ligand axes' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the central atom' },
      { id: '6S4', symbol: '6S4', operationType: 'Sn', count: 6, description: '90° improper rotations' },
      { id: '8S6', symbol: '8S6', operationType: 'Sn', count: 8, description: '60° improper rotations' },
      { id: '3sigma_h', symbol: '3σh', operationType: 'sigma_h', count: 3, description: 'Three Cartesian coordinate mirror planes (xy, xz, yz)' },
      { id: '6sigma_d', symbol: '6σd', operationType: 'sigma_d', count: 6, description: 'Six diagonal mirror planes bisecting the Cartesian axes' },
    ],
    irreps: [
      { symbol: 'A1g', id: 'A1g', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²+z²'], orbitalBasis: 's', irActive: false, ramanActive: true, description: 'Totally symmetric spherical representation.' },
      { symbol: 'Eg', id: 'Eg', dimension: 2, characters: [2, -1, 0, 0, 2, 2, 0, -1, 2, 0], linearFunctions: [], quadraticFunctions: ['(2z²-x²-y², x²-y²)'], orbitalBasis: '(dz², dx²-y²)', irActive: false, ramanActive: true, description: 'Doubly degenerate axial d-orbitals pointing directly at the 6 octahedral ligands.' },
      { symbol: 'T2g', id: 'T2g', dimension: 3, characters: [3, 0, 1, -1, -1, 3, -1, 0, -1, 1], linearFunctions: [], quadraticFunctions: ['(xy, xz, yz)'], orbitalBasis: '(dxy, dxz, dyz)', irActive: false, ramanActive: true, description: 'Triply degenerate inter-axial d-orbitals (π symmetry in octahedra).' },
      { symbol: 'T1u', id: 'T1u', dimension: 3, characters: [3, 0, -1, 1, -1, -3, -1, 0, 1, 1], linearFunctions: ['(x, y, z)'], quadraticFunctions: [], orbitalBasis: '(px, py, pz)', irActive: true, ramanActive: false, description: 'Triply degenerate ungerade p-orbitals along x, y, z axes.' },
      { symbol: 'T2u', id: 'T2u', dimension: 3, characters: [3, 0, 1, -1, -1, -3, 1, 0, 1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'Non-bonding LGOs', irActive: false, ramanActive: false, description: 'Triply degenerate ungerade ligand group representation.' },
    ],
    elements3D: [
      { id: 'c4_z', type: 'axis', label: 'C4 (z-axis)', schoenflies: 'C4', fold: 4, axis: [0, 1, 0], color: '#2563EB', description: '4-fold rotation axis along trans S–F bonds' },
      { id: 'c3_diag', type: 'axis', label: 'C3 (Face Diagonal)', schoenflies: 'C3', fold: 3, axis: [0.577, 0.577, 0.577], color: '#0284C7', description: '3-fold axis passing through opposite octahedral faces' },
      { id: 'sigma_h1', type: 'plane', label: 'σh (Equatorial Plane)', schoenflies: 'σh', axis: [0, 1, 0], color: '#D97706', description: 'Horizontal mirror plane containing four equatorial F atoms' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Central sulfur atom inversion center' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_6', label: 'Γ_σ (6 Octahedral S–F σ Bonds)', characters: [6, 0, 0, 2, 2, 0, 0, 0, 4, 2], description: 'Six octahedral σ bonds decompose into Γ_σ = A1g + Eg + T1u' },
    ],
  },

  Cinfv: {
    id: 'Cinfv',
    schoenflies: 'C∞v',
    name: 'Linear Heteronuclear',
    order: 6,
    family: 'Linear',
    description: 'Linear molecules without an inversion center: infinite rotation axis C∞ along the internuclear z-axis and an infinite number of vertical mirror planes σv containing the bond axis. Examples: CO, HF, HCl, HCN, N2O.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2Cinf', symbol: '2C∞(φ)', operationType: 'Cn', count: 2, description: 'Rotation by arbitrary angle φ around the internuclear z-axis' },
      { id: 'sigma_v', symbol: '∞σv', operationType: 'sigma_v', count: 3, description: 'Any vertical mirror plane containing the internuclear z-axis' },
    ],
    irreps: [
      { symbol: 'Σ+', id: 'SigmaPlus', dimension: 1, characters: [1, 1, 1], linearFunctions: ['z'], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, pz, dz² (σ orbitals)', irActive: true, ramanActive: true, description: 'Cylindrically symmetric σ orbitals (Λ = 0) symmetric under σv reflection.' },
      { symbol: 'Σ-', id: 'SigmaMinus', dimension: 1, characters: [1, 1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Cylindrically symmetric state antisymmetric under σv.' },
      { symbol: 'Π', id: 'Pi', dimension: 2, characters: [2, -1, 0], linearFunctions: ['(x, y)', '(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(px, py), (dxz, dyz) (π orbitals)', irActive: true, ramanActive: true, description: 'Doubly degenerate π orbitals (Λ = 1) with one nodal plane containing the bond axis.' },
      { symbol: 'Δ', id: 'Delta', dimension: 2, characters: [2, -1, 0], linearFunctions: [], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(dx²-y², dxy) (δ orbitals)', irActive: false, ramanActive: true, description: 'Doubly degenerate δ orbitals (Λ = 2) with two nodal planes.' },
    ],
    elements3D: [
      { id: 'cinf_z', type: 'axis', label: 'C∞ (Internuclear Axis)', schoenflies: 'C∞', fold: 8, axis: [1, 0, 0], color: '#2563EB', description: 'Continuous cylindrical rotation axis along the bond' },
      { id: 'sigma_v1', type: 'plane', label: 'σv(xy)', schoenflies: 'σv', axis: [0, 0, 1], color: '#0D9488', description: 'Vertical mirror plane containing the internuclear axis' },
      { id: 'sigma_v2', type: 'plane', label: 'σv(xz)', schoenflies: 'σv', axis: [0, 1, 0], color: '#D97706', description: 'Orthogonal vertical mirror plane containing the internuclear axis' },
    ],
    defaultReduciblePresets: [
      { id: 'valence_co', label: 'Γ_val (CO 8 Valence AOs)', characters: [8, 2, 4], description: 'C(2s,2p) + O(2s,2p) decompose into 4Σ+ + 2Π' },
    ],
  },

  Dinfh: {
    id: 'Dinfh',
    schoenflies: 'D∞h',
    name: 'Linear Centrosymmetric',
    order: 12,
    family: 'Linear',
    description: 'Linear molecules with a center of inversion i: C∞ axis, ∞ perpendicular C2\' axes, horizontal mirror plane σh, and ∞ vertical planes σv. Examples: CO2, C2H2 (Acetylene), N2, O2, H2, F2.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2Cinf', symbol: '2C∞(φ)', operationType: 'Cn', count: 2, description: 'Rotation by angle φ around the internuclear z-axis' },
      { id: 'sigma_v', symbol: '∞σv', operationType: 'sigma_v', count: 3, description: 'Vertical mirror planes containing the internuclear axis' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the molecular center' },
      { id: '2Sinf', symbol: '2S∞(φ)', operationType: 'Sn', count: 2, description: 'Improper rotation around the internuclear axis' },
      { id: 'C2p', symbol: '∞C2\'', operationType: 'Cn', count: 3, description: 'Perpendicular 180° axes passing through the center' },
    ],
    irreps: [
      { symbol: 'Σg+', id: 'SigmagPlus', dimension: 1, characters: [1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz² (σg)', irActive: false, ramanActive: true, description: 'Centrosymmetric gerade σ orbital.' },
      { symbol: 'Πg', id: 'Pig', dimension: 2, characters: [2, -1, 0, 2, -1, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz) (πg non-bonding / π*)', irActive: false, ramanActive: true, description: 'Doubly degenerate gerade π orbital (HOMO of CO2 / LUMO of C2H2).' },
      { symbol: 'Σu+', id: 'SigmauPlus', dimension: 1, characters: [1, 1, 1, -1, -1, -1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz (σu)', irActive: true, ramanActive: false, description: 'Antisymmetric ungerade σ orbital along the molecular axis.' },
      { symbol: 'Πu', id: 'Piu', dimension: 2, characters: [2, -1, 0, -2, 1, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py) (πu bonding)', irActive: true, ramanActive: false, description: 'Doubly degenerate ungerade π bonding orbital.' },
    ],
    elements3D: [
      { id: 'cinf_z', type: 'axis', label: 'C∞ (Molecular Axis)', schoenflies: 'C∞', fold: 8, axis: [1, 0, 0], color: '#2563EB', description: 'Continuous rotation axis along the linear backbone' },
      { id: 'c2_perp', type: 'axis', label: 'C2\' (Perpendicular)', schoenflies: 'C2\'', fold: 2, axis: [0, 1, 0], color: '#0284C7', description: 'Perpendicular 2-fold rotation axis' },
      { id: 'sigma_h', type: 'plane', label: 'σh (Central Plane)', schoenflies: 'σh', axis: [1, 0, 0], color: '#D97706', description: 'Horizontal mirror plane bisecting the linear molecule' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion' },
    ],
    defaultReduciblePresets: [
      { id: 'co2_sigma', label: 'Γ_σ (2 C=O σ Bonds)', characters: [2, 2, 2, 0, 0, 0], description: 'Two collinear C–O σ bonds decompose into Σg+ + Σu+' },
    ],
  },

  Cs: {
    id: 'Cs',
    schoenflies: 'Cs',
    name: 'Planar Reflection Symmetry',
    order: 2,
    family: 'Low Symmetry',
    description: 'Possesses only the identity operation E and a single mirror plane σh (usually the molecular plane). All orbitals classify into in-plane symmetric (A\') or out-of-plane antisymmetric π (A\'\').',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: 'sigma_h', symbol: 'σh(xy)', operationType: 'sigma_h', count: 1, description: 'Reflection across the molecular plane' },
    ],
    irreps: [
      { symbol: 'A\'', id: 'Ap', dimension: 1, characters: [1, 1], linearFunctions: ['x', 'y', 'Rz'], quadraticFunctions: ['x²', 'y²', 'z²', 'xy'], orbitalBasis: 's, px, py (σ framework)', irActive: true, ramanActive: true, description: 'Symmetric with respect to reflection in the molecular plane (σ orbitals).' },
      { symbol: 'A\'\'', id: 'App', dimension: 1, characters: [1, -1], linearFunctions: ['z', 'Rx', 'Ry'], quadraticFunctions: ['xz', 'yz'], orbitalBasis: 'pz (π framework)', irActive: true, ramanActive: true, description: 'Antisymmetric with respect to reflection in the molecular plane (π orbitals).' },
    ],
    elements3D: [
      { id: 'sigma_h', type: 'plane', label: 'σh (Molecular Plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#0D9488', description: 'Mirror plane containing the planar molecular skeleton' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_pi', label: 'Γ_val (4 In-plane A\' + 2 Out-of-plane A\'\')', characters: [6, 2], description: 'Decomposes into 4A\' (σ) + 2A\'\' (π)' },
    ],
  },
};

/**
 * Applies the Great Orthogonality Theorem to decompose a reducible representation
 * Γ_red into a direct sum of Irreducible Representations:
 *   a_i = (1 / h) * sum_k (N_k * χ_red(R_k) * χ_i(R_k))
 */
export function decomposeReducibleRepresentation(
  pointGroup: PointGroupData,
  reducibleChars: number[]
): {
  formula: string;
  results: DecompositionResult[];
  isValid: boolean;
} {
  const h = pointGroup.classes.reduce((acc, c) => acc + c.count, 0) || pointGroup.order;
  const results: DecompositionResult[] = [];
  const formulaParts: string[] = [];

  for (const irrep of pointGroup.irreps) {
    let sumContribution = 0;
    const terms: DecompositionResult['terms'] = [];

    for (let k = 0; k < pointGroup.classes.length; k++) {
      const cls = pointGroup.classes[k];
      const redChar = reducibleChars[k] ?? 0;
      const irrepChar = irrep.characters[k] ?? 0;
      const product = cls.count * redChar * irrepChar;
      sumContribution += product;
      terms.push({
        classSymbol: cls.symbol,
        count: cls.count,
        reducibleChar: redChar,
        irrepChar,
        product,
      });
    }

    const rawCoeff = sumContribution / h;
    const coefficient = Math.round(rawCoeff * 100) / 100;

    results.push({
      irrepSymbol: irrep.symbol,
      irrepId: irrep.id,
      coefficient,
      sumContribution,
      terms,
    });

    if (coefficient > 0.05) {
      const roundedInt = Math.abs(coefficient - Math.round(coefficient)) < 0.15
        ? Math.round(coefficient)
        : coefficient;
      formulaParts.push(roundedInt === 1 ? irrep.symbol : `${roundedInt}${irrep.symbol}`);
    }
  }

  return {
    formula: formulaParts.length > 0 ? formulaParts.join(' + ') : '0',
    results,
    isValid: formulaParts.length > 0,
  };
}
