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

  C4v: {
    id: 'C4v',
    schoenflies: 'C4v',
    name: 'Square Pyramidal',
    order: 8,
    family: 'Cnv (Pyramidal / Polar)',
    description: 'Features a 4-fold principal C4 axis, a collinear C2(=C4²) axis, two vertical mirror planes σv along the basal bonds, and two dihedral mirror planes σd bisecting the basal bonds. Classic examples: BrF5, IF5, XeOF4, SF5⁻.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C4', symbol: '2C4(z)', operationType: 'Cn', count: 2, description: '±90° rotations around the principal z-axis' },
      { id: 'C2', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation around the principal z-axis' },
      { id: '2sigma_v', symbol: '2σv', operationType: 'sigma_v', count: 2, description: 'Two vertical mirror planes containing the trans basal ligand pairs' },
      { id: '2sigma_d', symbol: '2σd', operationType: 'sigma_d', count: 2, description: 'Two dihedral mirror planes bisecting the basal bond angles' },
    ],
    irreps: [
      { symbol: 'A1', id: 'A1', dimension: 1, characters: [1, 1, 1, 1, 1], linearFunctions: ['z'], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, pz, dz²', irActive: true, ramanActive: true, description: 'Totally symmetric representation along the C4 axial bond and trans lone pair.' },
      { symbol: 'A2', id: 'A2', dimension: 1, characters: [1, 1, 1, -1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Symmetric under C4/C2 rotations, antisymmetric under σv and σd.' },
      { symbol: 'B1', id: 'B1', dimension: 1, characters: [1, -1, 1, 1, -1], linearFunctions: [], quadraticFunctions: ['x²-y²'], orbitalBasis: 'dx²-y²', irActive: false, ramanActive: true, description: 'Basal ligand σ-bonding symmetry matching lobes along the x and y axes.' },
      { symbol: 'B2', id: 'B2', dimension: 1, characters: [1, -1, 1, -1, 1], linearFunctions: [], quadraticFunctions: ['xy'], orbitalBasis: 'dxy', irActive: false, ramanActive: true, description: 'In-plane dxy representation bisecting the basal bonds.' },
      { symbol: 'E', id: 'E', dimension: 2, characters: [2, 0, -2, 0, 0], linearFunctions: ['(x, y)', '(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(px, py), (dxz, dyz)', irActive: true, ramanActive: true, description: 'Doubly degenerate representation spanning the basal equatorial plane.' },
    ],
    elements3D: [
      { id: 'c4_z', type: 'axis', label: 'C4 (z-axis)', schoenflies: 'C4', fold: 4, axis: [0, 1, 0], color: '#2563EB', description: 'Principal 4-fold rotation axis along the axial bond' },
      { id: 'sigma_v1', type: 'plane', label: 'σv(xz)', schoenflies: 'σv', axis: [0, 0, 1], color: '#0D9488', description: 'Vertical mirror plane containing trans basal bonds' },
      { id: 'sigma_v2', type: 'plane', label: 'σv(yz)', schoenflies: 'σv', axis: [1, 0, 0], color: '#0D9488', description: 'Orthogonal vertical mirror plane containing trans basal bonds' },
      { id: 'sigma_d1', type: 'plane', label: 'σd (Dihedral)', schoenflies: 'σd', axis: [0.707, 0, -0.707], color: '#D97706', description: 'Dihedral mirror plane bisecting basal bond angles' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_5', label: 'Γ_σ (5 Square Pyramidal Bonds)', characters: [5, 1, 1, 3, 1], description: '1 axial (A1) + 4 basal (A1 + B1 + E) bonds decompose into Γ_σ = 2A1 + B1 + E' },
    ],
  },

  D2d: {
    id: 'D2d',
    schoenflies: 'D2d',
    name: 'Tetragonal Disphenoidal / Bis-Diene Chelate',
    order: 8,
    family: 'Dnd (Dihedral Antiprismatic)',
    description: 'Contains an improper 4-fold S4 axis, a collinear C2(z) axis, two perpendicular C2\' axes, and two dihedral mirror planes σd bisecting the C2\' axes. Classic examples: Bis(1,5-cyclooctadiene)nickel(0) Ni(COD)2, Allene (H2C=C=CH2), Cyclooctatetraene (tub C8H8).',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2S4', symbol: '2S4(z)', operationType: 'Sn', count: 2, description: '±90° improper rotary-reflections around the principal S4 axis' },
      { id: 'C2z', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation around the principal S4/C2 axis' },
      { id: '2C2p', symbol: '2C2\'', operationType: 'Cn', count: 2, description: 'Two perpendicular 180° rotation axes orthogonal to S4' },
      { id: '2sigma_d', symbol: '2σd', operationType: 'sigma_d', count: 2, description: 'Two dihedral mirror planes containing the S4 axis and bisecting the 2C2\' axes' },
    ],
    irreps: [
      { symbol: 'A1', id: 'A1', dimension: 1, characters: [1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation invariant under S4, C2, and σd.' },
      { symbol: 'A2', id: 'A2', dimension: 1, characters: [1, 1, 1, -1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Symmetric under S4 and C2(z), antisymmetric under C2\' and σd.' },
      { symbol: 'B1', id: 'B1', dimension: 1, characters: [1, -1, 1, 1, -1], linearFunctions: [], quadraticFunctions: ['x²-y²'], orbitalBasis: 'dx²-y²', irActive: false, ramanActive: true, description: 'Antisymmetric under S4 and σd, symmetric under C2(z) and C2\'.' },
      { symbol: 'B2', id: 'B2', dimension: 1, characters: [1, -1, 1, -1, 1], linearFunctions: ['z'], quadraticFunctions: ['xy'], orbitalBasis: 'pz, dxy', irActive: true, ramanActive: true, description: 'Dipole-allowed representation along the S4/C2(z) axis and dxy backbonding manifold.' },
      { symbol: 'E', id: 'E', dimension: 2, characters: [2, 0, -2, 0, 0], linearFunctions: ['(x, y)', '(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(px, py), (dxz, dyz)', irActive: true, ramanActive: true, description: 'Doubly degenerate representation spanning the perpendicular (xz, yz) π-backbonding planes.' },
    ],
    elements3D: [
      { id: 's4_x', type: 'axis', label: 'S4 / C2 (Principal Axis)', schoenflies: 'S4', fold: 4, axis: [1, 0, 0], color: '#2563EB', description: 'Improper 4-fold S4 and collinear C2 axis passing through both chelating diene centroids and the metal center' },
      { id: 'c2_p1', type: 'axis', label: 'C2\' (Bisector 1)', schoenflies: 'C2\'', fold: 2, axis: [0, 0.707, 0.707], color: '#0284C7', description: 'Perpendicular 2-fold axis exchanging the two orthogonal diene ligands' },
      { id: 'sigma_d1', type: 'plane', label: 'σd(xy)', schoenflies: 'σd', axis: [0, 0, 1], color: '#0D9488', description: 'Dihedral mirror plane containing the S4 axis and bisecting the C=C bonds of one diene' },
      { id: 'sigma_d2', type: 'plane', label: 'σd(xz)', schoenflies: 'σd', axis: [0, 1, 0], color: '#D97706', description: 'Orthogonal dihedral mirror plane containing the S4 axis and bisecting the second diene' },
    ],
    defaultReduciblePresets: [
      { id: 'diene_pi_4', label: 'Γ_π (4 Olefin π-Donor Orbitals)', characters: [4, 0, 0, 0, 2], description: 'Four η²-alkene π-donor lobes in M(diene)2 decompose into Γ_π = A1 + B2 + E' },
      { id: 'gamma_d', label: 'Γ_d (5 Metal d Orbitals in D2d)', characters: [5, -1, 1, 1, 1], description: 'Five metal d orbitals split into A1(dz²) + B1(dx²-y²) + B2(dxy) + E(dxz, dyz)' },
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
      { id: 'valence_12', label: 'Γ_val (12 Ethylene Valence AOs)', characters: [12, 0, 0, 0, 0, 0, 4, 8], description: 'All 12 valence AOs of C2H4 (2×C 2s,2px,2py,2pz + 4×H 1s) decompose into 3Ag + B2g + 2B3g + 3B1u + 2B2u + B3u' },
      { id: 'pi_system', label: 'Γ_π (2 Carbon 2px π Orbitals)', characters: [2, -2, 0, 0, 0, 2, 0, -2], description: 'Two out-of-plane 2px orbitals form bonding B3u (π) and antibonding B2g (π*)' },
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

  D3: {
    id: 'D3',
    schoenflies: 'D3',
    name: 'Chiral Dihedral / Tris-Chelate Propeller',
    order: 6,
    family: 'Dn (Chiral Dihedral)',
    description: 'Features a principal 3-fold C3 axis and three perpendicular C2\' axes with no mirror planes or inversion center (chiral Λ/Δ propeller). Classic examples: [Co(en)3]³⁺, [Ru(bpy)3]²⁺, Fe(acac)3, [Fe(C2O4)3]³⁻.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C3', symbol: '2C3(z)', operationType: 'Cn', count: 2, description: '±120° rotations around the principal 3-fold propeller axis' },
      { id: '3C2p', symbol: '3C2\'', operationType: 'Cn', count: 3, description: 'Three 180° rotations bisecting each bidentate chelate ring' },
    ],
    irreps: [
      { symbol: 'A1', id: 'A1', dimension: 1, characters: [1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation along the C3 propeller axis.' },
      { symbol: 'A2', id: 'A2', dimension: 1, characters: [1, 1, -1], linearFunctions: ['z', 'Rz'], quadraticFunctions: [], orbitalBasis: 'pz', irActive: true, ramanActive: false, description: 'Symmetric under C3 rotation, antisymmetric under perpendicular C2\' rotations.' },
      { symbol: 'E', id: 'E', dimension: 2, characters: [2, -1, 0], linearFunctions: ['(x, y)', '(Rx, Ry)'], quadraticFunctions: ['(x²-y², xy)', '(xz, yz)'], orbitalBasis: '(px, py), (dxz, dyz), (dx²-y², dxy)', irActive: true, ramanActive: true, description: 'Doubly degenerate representation spanning the plane perpendicular to C3.' },
    ],
    elements3D: [
      { id: 'c3_diag', type: 'axis', label: 'C3 (Propeller Axis)', schoenflies: 'C3', fold: 3, axis: [0.577, 0.577, 0.577], color: '#2563EB', description: 'Principal 3-fold axis of the tris-chelate propeller' },
      { id: 'c2_1', type: 'axis', label: 'C2\' (Chelate Axis)', schoenflies: 'C2\'', fold: 2, axis: [0.707, 0.707, 0], color: '#0284C7', description: '2-fold rotation axis bisecting a bidentate chelate ligand' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_6', label: 'Γ_σ (6 Tris-Chelate M–L σ Bonds)', characters: [6, 0, 0], description: 'Six octahedral chelate σ donor bonds decompose into Γ_σ = A1 + A2 + 2E' },
      { id: 'gamma_d', label: 'Γ_d (5 Metal d Orbitals in D3)', characters: [5, -1, 1], description: 'Octahedral t2g splits into A1 + E while eg remains E (Γ_d = A1 + 2E)' },
    ],
  },

  D5d: {
    id: 'D5d',
    schoenflies: 'D5d',
    name: 'Pentagonal Antiprismatic / Staggered Metallocene',
    order: 20,
    family: 'Dnd (Dihedral Antiprismatic)',
    description: 'Contains a 5-fold principal C5 axis, five perpendicular C2\' axes, an inversion center i, an S10 improper axis, and five dihedral mirror planes σd. Classic examples: Staggered Ferrocene Fe(C5H5)2, Cobaltocene Co(C5H5)2, Ruthenocene Ru(C5H5)2.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C5', symbol: '2C5(z)', operationType: 'Cn', count: 2, description: '±72° rotations around the principal 5-fold sandwich axis' },
      { id: '2C5_2', symbol: '2C5²(z)', operationType: 'Cn', count: 2, description: '±144° rotations around the principal 5-fold axis' },
      { id: '5C2p', symbol: '5C2\'', operationType: 'Cn', count: 5, description: 'Five perpendicular 180° axes in the xy plane' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the central metal atom' },
      { id: '2S10_3', symbol: '2S10³', operationType: 'Sn', count: 2, description: '±108° improper rotations around the z-axis' },
      { id: '2S10', symbol: '2S10', operationType: 'Sn', count: 2, description: '±36° improper rotations around the z-axis' },
      { id: '5sigma_d', symbol: '5σd', operationType: 'sigma_d', count: 5, description: 'Five vertical dihedral mirror planes bisecting opposite C–C bonds' },
    ],
    irreps: [
      { symbol: 'A1g', id: 'A1g', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric metal dz² and Cp ring a1g π combination.' },
      { symbol: 'A2g', id: 'A2g', dimension: 1, characters: [1, 1, 1, -1, 1, 1, 1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'In-plane ring rotation.' },
      { symbol: 'E1g', id: 'E1g', dimension: 2, characters: [2, 0.618, -1.618, 0, 2, 0.618, -1.618, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz) — Strong M–Cp π Bonding', irActive: false, ramanActive: true, description: 'Primary covalent metal–cyclopentadienyl dxz/dyz π-bonding pair.' },
      { symbol: 'E2g', id: 'E2g', dimension: 2, characters: [2, -1.618, 0.618, 0, 2, -1.618, 0.618, 0], linearFunctions: [], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(dx²-y², dxy) — Weak δ Bonding', irActive: false, ramanActive: true, description: 'In-plane metal dx²-y²/dxy δ-symmetry frontier orbital pair.' },
      { symbol: 'A2u', id: 'A2u', dimension: 1, characters: [1, 1, 1, -1, -1, -1, -1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz', irActive: true, ramanActive: false, description: 'Antisymmetric out-of-plane metal pz + Cp ring combination.' },
      { symbol: 'E1u', id: 'E1u', dimension: 2, characters: [2, 0.618, -1.618, 0, -2, -0.618, 1.618, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py)', irActive: true, ramanActive: false, description: 'Ungerade metal px/py + Cp ring e1u combination.' },
    ],
    elements3D: [
      { id: 'c5_z', type: 'axis', label: 'C5 (z-axis)', schoenflies: 'C5', fold: 5, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 5-fold axis passing through both Cp ring centroids and the central metal' },
      { id: 'c2_p1', type: 'axis', label: 'C2\' (xy plane)', schoenflies: 'C2\'', fold: 2, axis: [1, 0, 0], color: '#0284C7', description: 'Perpendicular 2-fold axis through the central metal' },
      { id: 'sigma_d1', type: 'plane', label: 'σd (Dihedral)', schoenflies: 'σd', axis: [0, 1, 0], color: '#0D9488', description: 'Dihedral mirror plane containing the C5 axis' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Inversion center at the metal nucleus' },
    ],
    defaultReduciblePresets: [
      { id: 'cp_pi_10', label: 'Γ_π (10 Cp Carbon 2pz Orbitals)', characters: [10, 0, 0, 0, 0, 0, 0, 2], description: 'Ten Cp π orbitals decompose into A1g + A2u + E1g + E1u + E2g + E2u' },
      { id: 'gamma_d', label: 'Γ_d (5 Metal d Orbitals in D5d)', characters: [5, 0, 0, 1, 5, 0, 0, 1], description: 'Five metal d orbitals split into A1g(dz²) + E1g(dxz,dyz) + E2g(dx²-y²,dxy)' },
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
      { symbol: 'T1', id: 'T1', dimension: 3, characters: [3, 0, -1, 1, -1], linearFunctions: ['(Rx, Ry, Rz)'], quadraticFunctions: [], orbitalBasis: '(fx(x²-3y²), fy(3x²-y²), fz(x²-y²))', irActive: false, ramanActive: false, description: 'Triply degenerate rotational and cubic f-orbital representation.' },
      { symbol: 'T2', id: 'T2', dimension: 3, characters: [3, 0, -1, -1, 1], linearFunctions: ['(x, y, z)'], quadraticFunctions: ['(xy, xz, yz)'], orbitalBasis: '(px, py, pz), (dxy, dxz, dyz), (fxz², fyz², fz³)', irActive: true, ramanActive: true, description: 'Triply degenerate p, d, and axial f orbitals (sp³ / sd³ tetrahedral bonding).' },
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
    description: 'Full octahedral symmetry with inversion center i, three 4-fold C4 axes along the Cartesian x, y, z ligand axes, four 3-fold C3 axes through octahedral faces, and nine mirror planes (3 σh + 6 σd). Classic examples: SF6, [CrCl6]3-, [CeCl6]3-, Mo(CO)6.',
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
      { symbol: 'A2u', id: 'A2u', dimension: 1, characters: [1, 1, -1, -1, 1, -1, 1, -1, -1, 1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'fxyz', irActive: false, ramanActive: false, description: 'Cubic 8-lobed fxyz representation pointing into the eight octants between octahedral axes.' },
      { symbol: 'T1u', id: 'T1u', dimension: 3, characters: [3, 0, -1, 1, -1, -3, -1, 0, 1, 1], linearFunctions: ['(x, y, z)'], quadraticFunctions: [], orbitalBasis: '(px, py, pz), (fxz², fyz², fz³)', irActive: true, ramanActive: false, description: 'Triply degenerate ungerade p and axial f orbitals along x, y, z axes.' },
      { symbol: 'T2u', id: 'T2u', dimension: 3, characters: [3, 0, 1, -1, -1, -3, 1, 0, 1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '(fz(x²-y²), fx(x²-3y²), fy(3x²-y²))', irActive: false, ramanActive: false, description: 'Triply degenerate ungerade cubic f-orbital and ligand group representation.' },
    ],
    elements3D: [
      { id: 'c4_z', type: 'axis', label: 'C4 (z-axis)', schoenflies: 'C4', fold: 4, axis: [0, 1, 0], color: '#2563EB', description: '4-fold rotation axis along trans S–F bonds' },
      { id: 'c3_diag', type: 'axis', label: 'C3 (Face Diagonal)', schoenflies: 'C3', fold: 3, axis: [0.577, 0.577, 0.577], color: '#0284C7', description: '3-fold axis passing through opposite octahedral faces' },
      { id: 'sigma_h1', type: 'plane', label: 'σh (Equatorial Plane)', schoenflies: 'σh', axis: [0, 1, 0], color: '#D97706', description: 'Horizontal mirror plane containing four equatorial F atoms' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Central sulfur atom inversion center' },
    ],
    defaultReduciblePresets: [
      { id: 'sigma_6', label: 'Γ_σ (6 Octahedral M–L σ Bonds)', characters: [6, 0, 0, 2, 2, 0, 0, 0, 4, 2], description: 'Six octahedral σ bonds decompose into Γ_σ = A1g + Eg + T1u' },
      { id: 'gamma_d', label: 'Γ_d (5 Metal d Orbitals)', characters: [5, -1, 1, -1, 1, 5, -1, -1, 1, 1], description: 'Five central metal d orbitals split into Γ_d = Eg + T2g' },
      { id: 'gamma_f', label: 'Γ_f (7 Lanthanide 4f Orbitals)', characters: [7, 1, -1, -1, -1, -7, 1, -1, 1, 1], description: 'Seven central lanthanide 4f orbitals split into Γ_f = A2u + T1u + T2u' },
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

  C1: {
    id: 'C1',
    schoenflies: 'C1',
    name: 'Trivial / Chiral Asymmetric',
    order: 1,
    family: 'Low Symmetry',
    description: 'Contains only the identity operation E. Characteristic of chiral asymmetric molecules (e.g., CHFClBr, ATP, amino acids). All orbitals transform as the totally symmetric A representation.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
    ],
    irreps: [
      { symbol: 'A', id: 'A', dimension: 1, characters: [1], linearFunctions: ['x', 'y', 'z', 'Rx', 'Ry', 'Rz'], quadraticFunctions: ['x²', 'y²', 'z²', 'xy', 'xz', 'yz'], orbitalBasis: 'All s, p, d, f orbitals', irActive: true, ramanActive: true, description: 'Sole totally symmetric irreducible representation of the trivial group.' },
    ],
    elements3D: [],
    defaultReduciblePresets: [
      { id: 'c1_val', label: 'Γ_val (4 Valence AOs)', characters: [4], description: 'Trivially decomposes into 4A' },
    ],
  },

  Ci: {
    id: 'Ci',
    schoenflies: 'Ci',
    name: 'Centrosymmetric Inversion',
    order: 2,
    family: 'Low Symmetry',
    description: 'Contains only the identity E and a center of inversion i. Orbitals partition strictly by parity into gerade (Ag) and ungerade (Au). Examples: meso-tartaric acid, anti-1,2-dichloro-1,2-difluoroethane.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the molecular center' },
    ],
    irreps: [
      { symbol: 'Ag', id: 'Ag', dimension: 1, characters: [1, 1], linearFunctions: ['Rx', 'Ry', 'Rz'], quadraticFunctions: ['x²', 'y²', 'z²', 'xy', 'xz', 'yz'], orbitalBasis: 's, d-orbitals (gerade)', irActive: false, ramanActive: true, description: 'Symmetric (gerade) under inversion; Raman active, IR forbidden (mutual exclusion).' },
      { symbol: 'Au', id: 'Au', dimension: 1, characters: [1, -1], linearFunctions: ['x', 'y', 'z'], quadraticFunctions: [], orbitalBasis: 'px, py, pz, f-orbitals (ungerade)', irActive: true, ramanActive: false, description: 'Antisymmetric (ungerade) under inversion; IR active, Raman forbidden.' },
    ],
    elements3D: [
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion' },
    ],
    defaultReduciblePresets: [
      { id: 'ci_pair', label: 'Γ_pair (2 Inversion-Related AOs)', characters: [2, 0], description: 'Decomposes into Ag + Au' },
    ],
  },

  C2: {
    id: 'C2',
    schoenflies: 'C2',
    name: 'Chiral 2-Fold Rotational',
    order: 2,
    family: 'Cn (Chiral Axial)',
    description: 'Contains a single proper 2-fold rotation axis C2 without mirror planes or inversion. Chiral point group of H2O2 (gauche), hydrazine, and skewed biaryls.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: 'C2', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation around the principal C2 axis' },
    ],
    irreps: [
      { symbol: 'A', id: 'A', dimension: 1, characters: [1, 1], linearFunctions: ['z', 'Rz'], quadraticFunctions: ['x²', 'y²', 'z²', 'xy'], orbitalBasis: 's, pz, dz², dx²-y², dxy', irActive: true, ramanActive: true, description: 'Symmetric under C2(z) rotation.' },
      { symbol: 'B', id: 'B', dimension: 1, characters: [1, -1], linearFunctions: ['x', 'y', 'Rx', 'Ry'], quadraticFunctions: ['xz', 'yz'], orbitalBasis: 'px, py, dxz, dyz', irActive: true, ramanActive: true, description: 'Antisymmetric under C2(z) rotation.' },
    ],
    elements3D: [
      { id: 'c2_z', type: 'axis', label: 'C2 (z-axis)', schoenflies: 'C2', fold: 2, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 180° rotational symmetry axis' },
    ],
    defaultReduciblePresets: [
      { id: 'c2_bonds', label: 'Γ_σ (2 Equivalent O–H Bonds in H2O2)', characters: [2, 0], description: 'Decomposes into A + B' },
    ],
  },

  C2h: {
    id: 'C2h',
    schoenflies: 'C2h',
    name: 'Centrosymmetric Planar 2-Fold',
    order: 4,
    family: 'Cnh (Centrosymmetric Planar)',
    description: 'Features a C2 axis perpendicular to a horizontal mirror plane σh and an inversion center i. Classic examples: trans-1,3-butadiene, trans-1,2-dichloroethene, trans-N2F2.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: 'C2', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation perpendicular to molecular plane' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion through the molecular center' },
      { id: 'sigma_h', symbol: 'σh(xy)', operationType: 'sigma_h', count: 1, description: 'Horizontal mirror plane containing the molecule' },
    ],
    irreps: [
      { symbol: 'Ag', id: 'Ag', dimension: 1, characters: [1, 1, 1, 1], linearFunctions: ['Rz'], quadraticFunctions: ['x²', 'y²', 'z²', 'xy'], orbitalBasis: 's, dz², dx²-y², dxy', irActive: false, ramanActive: true, description: 'Totally symmetric in-plane gerade orbital.' },
      { symbol: 'Bg', id: 'Bg', dimension: 1, characters: [1, -1, 1, -1], linearFunctions: ['Rx', 'Ry'], quadraticFunctions: ['xz', 'yz'], orbitalBasis: 'dxz, dyz (π gerade)', irActive: false, ramanActive: true, description: 'Out-of-plane π gerade orbital (HOMO of trans-butadiene).' },
      { symbol: 'Au', id: 'Au', dimension: 1, characters: [1, 1, -1, -1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz (π ungerade)', irActive: true, ramanActive: false, description: 'Out-of-plane π ungerade orbital (lowest π bonding MO of trans-butadiene).' },
      { symbol: 'Bu', id: 'Bu', dimension: 1, characters: [1, -1, -1, 1], linearFunctions: ['x', 'y'], quadraticFunctions: [], orbitalBasis: 'px, py (σ ungerade)', irActive: true, ramanActive: false, description: 'In-plane σ ungerade orbital.' },
    ],
    elements3D: [
      { id: 'c2_z', type: 'axis', label: 'C2 (z-axis)', schoenflies: 'C2', fold: 2, axis: [0, 0, 1], color: '#2563EB', description: 'Principal C2 axis perpendicular to the molecular plane' },
      { id: 'sigma_h', type: 'plane', label: 'σh (Molecular Plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#0D9488', description: 'Horizontal mirror plane' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion' },
    ],
    defaultReduciblePresets: [
      { id: 'c2h_pi', label: 'Γ_π (trans-Butadiene 4 pz AOs)', characters: [4, 0, 0, -4], description: 'Decomposes into 2Au + 2Bg' },
    ],
  },

  C3h: {
    id: 'C3h',
    schoenflies: 'C3h',
    name: 'Planar 3-Fold Horizontal Mirror',
    order: 6,
    family: 'Cnh',
    description: 'Contains a 3-fold rotation axis C3, horizontal mirror plane σh, and improper S3 axis, without vertical mirror planes. Classic example: boric acid B(OH)3.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C3', symbol: '2C3', operationType: 'Cn', count: 2, description: '±120° rotations around principal axis' },
      { id: 'sigma_h', symbol: 'σh', operationType: 'sigma_h', count: 1, description: 'Horizontal mirror plane' },
      { id: '2S3', symbol: '2S3', operationType: 'Sn', count: 2, description: '±120° improper rotations' },
    ],
    irreps: [
      { symbol: 'A\'', id: 'Ap', dimension: 1, characters: [1, 1, 1, 1], linearFunctions: ['Rz'], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric in-plane representation.' },
      { symbol: 'E\'', id: 'Ep', dimension: 2, characters: [2, -1, 2, -1], linearFunctions: ['(x, y)'], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(px, py), (dx²-y², dxy)', irActive: true, ramanActive: true, description: 'Doubly degenerate in-plane representation.' },
      { symbol: 'A\'\'', id: 'App', dimension: 1, characters: [1, 1, -1, -1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz', irActive: true, ramanActive: false, description: 'Out-of-plane totally symmetric with respect to C3.' },
      { symbol: 'E\'\'', id: 'Epp', dimension: 2, characters: [2, -1, -2, 1], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz)', irActive: false, ramanActive: true, description: 'Doubly degenerate out-of-plane representation.' },
    ],
    elements3D: [
      { id: 'c3_z', type: 'axis', label: 'C3 (z-axis)', schoenflies: 'C3', fold: 3, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 3-fold rotation axis' },
      { id: 'sigma_h', type: 'plane', label: 'σh (xy-plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#0D9488', description: 'Horizontal mirror plane' },
    ],
    defaultReduciblePresets: [
      { id: 'c3h_sigma', label: 'Γ_σ (3 B–O Bonds)', characters: [3, 0, 3, 0], description: 'Decomposes into A\' + E\'' },
    ],
  },

  D2: {
    id: 'D2',
    schoenflies: 'D2',
    name: 'Chiral Orthorhombic Dihedral',
    order: 4,
    family: 'Dn (Chiral Dihedral)',
    description: 'Contains three mutually orthogonal C2 rotation axes (x, y, z) and no mirror planes or inversion center. Classic example: twisted biphenyl or twisted ethylene (0° < θ < 90°).',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: 'C2z', symbol: 'C2(z)', operationType: 'Cn', count: 1, description: '180° rotation about z' },
      { id: 'C2y', symbol: 'C2(y)', operationType: 'Cn', count: 1, description: '180° rotation about y' },
      { id: 'C2x', symbol: 'C2(x)', operationType: 'Cn', count: 1, description: '180° rotation about x' },
    ],
    irreps: [
      { symbol: 'A', id: 'A', dimension: 1, characters: [1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²', 'y²', 'z²'], orbitalBasis: 's, dz², dx²-y²', irActive: false, ramanActive: true, description: 'Totally symmetric under all three C2 axes.' },
      { symbol: 'B1', id: 'B1', dimension: 1, characters: [1, 1, -1, -1], linearFunctions: ['z', 'Rz'], quadraticFunctions: ['xy'], orbitalBasis: 'pz, dxy', irActive: true, ramanActive: true, description: 'Symmetric under C2(z), antisymmetric under C2(x) and C2(y).' },
      { symbol: 'B2', id: 'B2', dimension: 1, characters: [1, -1, 1, -1], linearFunctions: ['y', 'Ry'], quadraticFunctions: ['xz'], orbitalBasis: 'py, dxz', irActive: true, ramanActive: true, description: 'Symmetric under C2(y), antisymmetric under C2(z) and C2(x).' },
      { symbol: 'B3', id: 'B3', dimension: 1, characters: [1, -1, -1, 1], linearFunctions: ['x', 'Rx'], quadraticFunctions: ['yz'], orbitalBasis: 'px, dyz', irActive: true, ramanActive: true, description: 'Symmetric under C2(x), antisymmetric under C2(z) and C2(y).' },
    ],
    elements3D: [
      { id: 'c2_z', type: 'axis', label: 'C2(z)', schoenflies: 'C2', fold: 2, axis: [0, 0, 1], color: '#2563EB', description: '2-fold rotation axis along z' },
      { id: 'c2_y', type: 'axis', label: 'C2(y)', schoenflies: 'C2', fold: 2, axis: [0, 1, 0], color: '#0D9488', description: '2-fold rotation axis along y' },
      { id: 'c2_x', type: 'axis', label: 'C2(x)', schoenflies: 'C2', fold: 2, axis: [1, 0, 0], color: '#D97706', description: '2-fold rotation axis along x' },
    ],
    defaultReduciblePresets: [
      { id: 'd2_four', label: 'Γ_4 (4 Equivalent Substituents)', characters: [4, 0, 0, 0], description: 'Decomposes into A + B1 + B2 + B3' },
    ],
  },

  D3d: {
    id: 'D3d',
    schoenflies: 'D3d',
    name: 'Trigonal Antiprismatic / Staggered',
    order: 12,
    family: 'Dnd (Antiprismatic)',
    description: 'Features a C3 principal axis, 3 perpendicular C2 axes, an inversion center i, S6 improper axis, and 3 dihedral mirror planes σd. Classic examples: staggered ethane (C2H6), chair cyclohexane (C6H12).',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C3', symbol: '2C3', operationType: 'Cn', count: 2, description: '±120° rotations around C3' },
      { id: '3C2', symbol: '3C2\'', operationType: 'Cn', count: 3, description: 'Three perpendicular C2 axes bisecting staggered bonds' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion center' },
      { id: '2S6', symbol: '2S6', operationType: 'Sn', count: 2, description: '±60° improper rotations' },
      { id: '3sigma_d', symbol: '3σd', operationType: 'sigma_d', count: 3, description: 'Three dihedral mirror planes' },
    ],
    irreps: [
      { symbol: 'A1g', id: 'A1g', dimension: 1, characters: [1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric gerade representation.' },
      { symbol: 'A2g', id: 'A2g', dimension: 1, characters: [1, 1, -1, 1, 1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: 'Rz', irActive: false, ramanActive: false, description: 'Gerade representation antisymmetric under C2\' and σd.' },
      { symbol: 'Eg', id: 'Eg', dimension: 2, characters: [2, -1, 0, 2, -1, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(x²-y², xy)', '(xz, yz)'], orbitalBasis: '(dx²-y², dxy), (dxz, dyz)', irActive: false, ramanActive: true, description: 'Doubly degenerate gerade representation (HOMO of staggered ethane).' },
      { symbol: 'A1u', id: 'A1u', dimension: 1, characters: [1, 1, 1, -1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'f(x³-3xy²)', irActive: false, ramanActive: false, description: 'Ungerade representation symmetric under C3 and C2\'.' },
      { symbol: 'A2u', id: 'A2u', dimension: 1, characters: [1, 1, -1, -1, -1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz, fz³', irActive: true, ramanActive: false, description: 'Axial ungerade representation.' },
      { symbol: 'Eu', id: 'Eu', dimension: 2, characters: [2, -1, 0, -2, 1, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py)', irActive: true, ramanActive: false, description: 'Doubly degenerate ungerade representation.' },
    ],
    elements3D: [
      { id: 'c3_z', type: 'axis', label: 'C3 / S6 (z-axis)', schoenflies: 'C3', fold: 3, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 3-fold rotation & S6 improper axis' },
      { id: 'c2_1', type: 'axis', label: 'C2\' (Perpendicular)', schoenflies: 'C2\'', fold: 2, axis: [1, 0, 0], color: '#0284C7', description: 'Perpendicular C2 axis' },
      { id: 'sigma_d', type: 'plane', label: 'σd (Dihedral Plane)', schoenflies: 'σd', axis: [0, 1, 0], color: '#0D9488', description: 'Dihedral mirror plane' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion' },
    ],
    defaultReduciblePresets: [
      { id: 'ethane_ch', label: 'Γ_CH (6 Staggered C–H Bonds)', characters: [6, 0, 0, 0, 0, 2], description: 'Decomposes into A1g + Eg + A2u + Eu' },
    ],
  },

  D4d: {
    id: 'D4d',
    schoenflies: 'D4d',
    name: 'Square Antiprismatic',
    order: 16,
    family: 'Dnd (Antiprismatic)',
    description: 'Features a C4 principal axis, S8 improper axis, 4 perpendicular C2\' axes, and 4 dihedral mirror planes σd. Classic examples: S8 crown ring, Mn2(CO)10 (staggered), [TaF8]3-.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2S8', symbol: '2S8', operationType: 'Sn', count: 2, description: '±45° improper rotations' },
      { id: '2C4', symbol: '2C4', operationType: 'Cn', count: 2, description: '±90° proper rotations' },
      { id: '2S83', symbol: '2S8³', operationType: 'Sn', count: 2, description: '±135° improper rotations' },
      { id: 'C2', symbol: 'C2', operationType: 'Cn', count: 1, description: '180° rotation about principal axis' },
      { id: '4C2p', symbol: '4C2\'', operationType: 'Cn', count: 4, description: 'Four perpendicular 2-fold axes' },
      { id: '4sigma_d', symbol: '4σd', operationType: 'sigma_d', count: 4, description: 'Four dihedral mirror planes' },
    ],
    irreps: [
      { symbol: 'A1', id: 'A1', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation.' },
      { symbol: 'A2', id: 'A2', dimension: 1, characters: [1, 1, 1, 1, 1, -1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: 'Rz', irActive: false, ramanActive: false, description: 'Antisymmetric under C2\' and σd.' },
      { symbol: 'B1', id: 'B1', dimension: 1, characters: [1, -1, 1, -1, 1, 1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Antisymmetric under S8 and σd.' },
      { symbol: 'B2', id: 'B2', dimension: 1, characters: [1, -1, 1, -1, 1, -1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz', irActive: true, ramanActive: false, description: 'Axial polar representation along z.' },
      { symbol: 'E1', id: 'E1', dimension: 2, characters: [2, 1.4142, 0, -1.4142, -2, 0, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py)', irActive: true, ramanActive: false, description: 'Doubly degenerate IR-active representation.' },
      { symbol: 'E2', id: 'E2', dimension: 2, characters: [2, 0, -2, 0, 2, 0, 0], linearFunctions: [], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(dx²-y², dxy)', irActive: false, ramanActive: true, description: 'Doubly degenerate equatorial d-orbital representation.' },
      { symbol: 'E3', id: 'E3', dimension: 2, characters: [2, -1.4142, 0, 1.4142, -2, 0, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz)', irActive: false, ramanActive: true, description: 'Doubly degenerate π-type d-orbital representation.' },
    ],
    elements3D: [
      { id: 'c4_z', type: 'axis', label: 'C4 / S8 (z-axis)', schoenflies: 'C4', fold: 4, axis: [0, 0, 1], color: '#2563EB', description: 'Principal C4 and S8 axis' },
      { id: 'sigma_d', type: 'plane', label: 'σd (Dihedral Plane)', schoenflies: 'σd', axis: [1, 0, 0], color: '#0D9488', description: 'Dihedral mirror plane' },
    ],
    defaultReduciblePresets: [
      { id: 'd4d_s8', label: 'Γ_8 (8 Vertices of Square Antiprism / S8)', characters: [8, 0, 0, 0, 0, 0, 0], description: 'Decomposes into A1 + B2 + E1 + E2 + E3' },
    ],
  },

  D5h: {
    id: 'D5h',
    schoenflies: 'D5h',
    name: 'Pentagonal Bipyramidal / Eclipsed Sandwich',
    order: 20,
    family: 'Dnh (Prismatic / Sandwich)',
    description: 'Features a 5-fold principal axis C5, 5 perpendicular C2 axes, horizontal mirror plane σh, and 5 vertical mirror planes σv. Classic examples: eclipsed ferrocene Fe(Cp)2, cyclopentadienyl anion Cp-, IF7.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2C5', symbol: '2C5', operationType: 'Cn', count: 2, description: '±72° rotations' },
      { id: '2C52', symbol: '2C5²', operationType: 'Cn', count: 2, description: '±144° rotations' },
      { id: '5C2', symbol: '5C2\'', operationType: 'Cn', count: 5, description: 'Five perpendicular C2 axes' },
      { id: 'sigma_h', symbol: 'σh', operationType: 'sigma_h', count: 1, description: 'Horizontal mirror plane' },
      { id: '2S5', symbol: '2S5', operationType: 'Sn', count: 2, description: '±72° improper rotations' },
      { id: '2S53', symbol: '2S5³', operationType: 'Sn', count: 2, description: '±216° improper rotations' },
      { id: '5sigma_v', symbol: '5σv', operationType: 'sigma_v', count: 5, description: 'Five vertical mirror planes' },
    ],
    irreps: [
      { symbol: 'A1\'', id: 'A1p', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation.' },
      { symbol: 'A2\'', id: 'A2p', dimension: 1, characters: [1, 1, 1, -1, 1, 1, 1, -1], linearFunctions: ['Rz'], quadraticFunctions: [], orbitalBasis: 'Rz', irActive: false, ramanActive: false, description: 'In-plane rotational representation.' },
      { symbol: 'E1\'', id: 'E1p', dimension: 2, characters: [2, 0.618, -1.618, 0, 2, 0.618, -1.618, 0], linearFunctions: ['(x, y)'], quadraticFunctions: [], orbitalBasis: '(px, py)', irActive: true, ramanActive: false, description: 'Doubly degenerate in-plane IR-active representation.' },
      { symbol: 'E2\'', id: 'E2p', dimension: 2, characters: [2, -1.618, 0.618, 0, 2, -1.618, 0.618, 0], linearFunctions: [], quadraticFunctions: ['(x²-y², xy)'], orbitalBasis: '(dx²-y², dxy)', irActive: false, ramanActive: true, description: 'Doubly degenerate δ-type equatorial representation.' },
      { symbol: 'A1\'\'', id: 'A1pp', dimension: 1, characters: [1, 1, 1, 1, -1, -1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Antisymmetric under σh and σv.' },
      { symbol: 'A2\'\'', id: 'A2pp', dimension: 1, characters: [1, 1, 1, -1, -1, -1, -1, 1], linearFunctions: ['z'], quadraticFunctions: [], orbitalBasis: 'pz (lowest Cp π)', irActive: true, ramanActive: false, description: 'Out-of-plane totally symmetric π orbital.' },
      { symbol: 'E1\'\'', id: 'E1pp', dimension: 2, characters: [2, 0.618, -1.618, 0, -2, -0.618, 1.618, 0], linearFunctions: ['(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(dxz, dyz) (Cp π HOMO)', irActive: false, ramanActive: true, description: 'Doubly degenerate out-of-plane π bonding orbital.' },
      { symbol: 'E2\'\'', id: 'E2pp', dimension: 2, characters: [2, -1.618, 0.618, 0, -2, 1.618, -0.618, 0], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'Cp π* LUMO', irActive: false, ramanActive: false, description: 'Doubly degenerate out-of-plane π* antibonding orbital.' },
    ],
    elements3D: [
      { id: 'c5_z', type: 'axis', label: 'C5 (z-axis)', schoenflies: 'C5', fold: 5, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 5-fold rotation axis' },
      { id: 'sigma_h', type: 'plane', label: 'σh (Equatorial Plane)', schoenflies: 'σh', axis: [0, 0, 1], color: '#D97706', description: 'Horizontal mirror plane' },
      { id: 'sigma_v', type: 'plane', label: 'σv (Vertical Plane)', schoenflies: 'σv', axis: [1, 0, 0], color: '#0D9488', description: 'Vertical mirror plane' },
    ],
    defaultReduciblePresets: [
      { id: 'cp_pi', label: 'Γ_π (Cyclopentadienyl 5 pz AOs)', characters: [5, 0, 0, -1, -5, 0, 0, 1], description: 'Decomposes into A2\'\' + E1\'\' + E2\'\'' },
    ],
  },

  S4: {
    id: 'S4',
    schoenflies: 'S4',
    name: '4-Fold Improper Rotational',
    order: 4,
    family: 'Sn (Improper Axial)',
    description: 'Contains an improper 4-fold rotation-reflection axis S4 (which generates C2 = S4²), without mirror planes or inversion. Achiral despite lacking σ and i. Example: 1,3,5,7-tetrafluorocyclooctatetraene.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '2S4', symbol: '2S4', operationType: 'Sn', count: 2, description: '±90° improper rotations' },
      { id: 'C2', symbol: 'C2', operationType: 'Cn', count: 1, description: '180° rotation (= S4²)' },
    ],
    irreps: [
      { symbol: 'A', id: 'A', dimension: 1, characters: [1, 1, 1], linearFunctions: ['Rz'], quadraticFunctions: ['x²+y²', 'z²'], orbitalBasis: 's, dz²', irActive: false, ramanActive: true, description: 'Totally symmetric representation.' },
      { symbol: 'B', id: 'B', dimension: 1, characters: [1, -1, 1], linearFunctions: ['z'], quadraticFunctions: ['x²-y²', 'xy'], orbitalBasis: 'pz, dx²-y², dxy', irActive: true, ramanActive: true, description: 'Antisymmetric under S4, symmetric under C2.' },
      { symbol: 'E', id: 'E', dimension: 2, characters: [2, 0, -2], linearFunctions: ['(x, y)', '(Rx, Ry)'], quadraticFunctions: ['(xz, yz)'], orbitalBasis: '(px, py), (dxz, dyz)', irActive: true, ramanActive: true, description: 'Doubly degenerate representation.' },
    ],
    elements3D: [
      { id: 's4_z', type: 'axis', label: 'S4 / C2 (z-axis)', schoenflies: 'S4', fold: 4, axis: [0, 0, 1], color: '#2563EB', description: 'Principal 4-fold improper rotation axis' },
    ],
    defaultReduciblePresets: [
      { id: 's4_four', label: 'Γ_4 (4 S4-Related Substituents)', characters: [4, 0, 0], description: 'Decomposes into A + B + E' },
    ],
  },

  Th: {
    id: 'Th',
    schoenflies: 'Th',
    name: 'Pyritohedral (Centrosymmetric Tetrahedral)',
    order: 24,
    family: 'Cubic / Polyhedral',
    description: 'Combines tetrahedral rotational symmetry T with an inversion center i and three orthogonal mirror planes σh. Classic examples: [Fe(H2O)6]2+ (oriented ligands), [Co(NO2)6]3-.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '8C3', symbol: '8C3', operationType: 'Cn', count: 8, description: 'Eight 3-fold body-diagonal rotations' },
      { id: '3C2', symbol: '3C2', operationType: 'Cn', count: 3, description: 'Three orthogonal C2 axes' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion center' },
      { id: '8S6', symbol: '8S6', operationType: 'Sn', count: 8, description: 'Eight 6-fold improper rotations' },
      { id: '3sigma_h', symbol: '3σh', operationType: 'sigma_h', count: 3, description: 'Three orthogonal Cartesian mirror planes' },
    ],
    irreps: [
      { symbol: 'Ag', id: 'Ag', dimension: 1, characters: [1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²+z²'], orbitalBasis: 's', irActive: false, ramanActive: true, description: 'Totally symmetric gerade representation.' },
      { symbol: 'Eg', id: 'Eg', dimension: 2, characters: [2, -1, 2, 2, -1, 2], linearFunctions: [], quadraticFunctions: ['(2z²-x²-y², x²-y²)'], orbitalBasis: '(dz², dx²-y²)', irActive: false, ramanActive: true, description: 'Doubly degenerate gerade eg d-orbitals.' },
      { symbol: 'Tg', id: 'Tg', dimension: 3, characters: [3, 0, -1, 3, 0, -1], linearFunctions: ['(Rx, Ry, Rz)'], quadraticFunctions: ['(xy, xz, yz)'], orbitalBasis: '(dxy, dxz, dyz)', irActive: false, ramanActive: true, description: 'Triply degenerate gerade t2g-like d-orbitals.' },
      { symbol: 'Au', id: 'Au', dimension: 1, characters: [1, 1, 1, -1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'fxyz', irActive: false, ramanActive: false, description: 'Singly degenerate ungerade representation.' },
      { symbol: 'Eu', id: 'Eu', dimension: 2, characters: [2, -1, 2, -2, 1, -2], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Doubly degenerate ungerade representation.' },
      { symbol: 'Tu', id: 'Tu', dimension: 3, characters: [3, 0, -1, -3, 0, 1], linearFunctions: ['(x, y, z)'], quadraticFunctions: [], orbitalBasis: '(px, py, pz)', irActive: true, ramanActive: false, description: 'Triply degenerate IR-active ungerade representation.' },
    ],
    elements3D: [
      { id: 'c2_z', type: 'axis', label: 'C2(z)', schoenflies: 'C2', fold: 2, axis: [0, 0, 1], color: '#2563EB', description: 'Cartesian C2 axis' },
      { id: 'c3_1', type: 'axis', label: 'C3 (1,1,1)', schoenflies: 'C3', fold: 3, axis: [0.577, 0.577, 0.577], color: '#0284C7', description: 'Body-diagonal C3 axis' },
      { id: 'sigma_h', type: 'plane', label: 'σh(xy)', schoenflies: 'σh', axis: [0, 0, 1], color: '#0D9488', description: 'Cartesian mirror plane' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion' },
    ],
    defaultReduciblePresets: [
      { id: 'th_6', label: 'Γ_σ (6 Octahedral Metal–Ligand Bonds)', characters: [6, 0, 2, 0, 0, 4], description: 'Decomposes into Ag + Eg + Tu' },
    ],
  },

  Ih: {
    id: 'Ih',
    schoenflies: 'Ih',
    name: 'Icosahedral (Centrosymmetric)',
    order: 120,
    family: 'Cubic / Polyhedral',
    description: 'Highest finite molecular point group (order h = 120): six C5 axes, ten C3 axes, fifteen C2 axes, inversion i, and fifteen mirror planes. Classic examples: C60 buckminsterfullerene, [B12H12]2-.',
    classes: [
      { id: 'E', symbol: 'E', operationType: 'E', count: 1, description: 'Identity operation' },
      { id: '12C5', symbol: '12C5', operationType: 'Cn', count: 12, description: 'Twelve ±72° rotations around 6 C5 axes' },
      { id: '12C52', symbol: '12C5²', operationType: 'Cn', count: 12, description: 'Twelve ±144° rotations around 6 C5 axes' },
      { id: '20C3', symbol: '20C3', operationType: 'Cn', count: 20, description: 'Twenty ±120° rotations' },
      { id: '15C2', symbol: '15C2', operationType: 'Cn', count: 15, description: 'Fifteen 180° rotations' },
      { id: 'i', symbol: 'i', operationType: 'i', count: 1, description: 'Inversion center' },
      { id: '12S10', symbol: '12S10', operationType: 'Sn', count: 12, description: 'Twelve improper S10 rotations' },
      { id: '12S103', symbol: '12S10³', operationType: 'Sn', count: 12, description: 'Twelve improper S10³ rotations' },
      { id: '20S6', symbol: '20S6', operationType: 'Sn', count: 20, description: 'Twenty improper S6 rotations' },
      { id: '15sigma', symbol: '15σ', operationType: 'sigma_v', count: 15, description: 'Fifteen mirror planes' },
    ],
    irreps: [
      { symbol: 'Ag', id: 'Ag', dimension: 1, characters: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], linearFunctions: [], quadraticFunctions: ['x²+y²+z²'], orbitalBasis: 's', irActive: false, ramanActive: true, description: 'Totally symmetric icosahedral representation.' },
      { symbol: 'T1g', id: 'T1g', dimension: 3, characters: [3, 1.618, -0.618, 0, -1, 3, -0.618, 1.618, 0, -1], linearFunctions: ['(Rx, Ry, Rz)'], quadraticFunctions: [], orbitalBasis: 'Rotations', irActive: false, ramanActive: false, description: 'Triply degenerate gerade representation.' },
      { symbol: 'T2g', id: 'T2g', dimension: 3, characters: [3, -0.618, 1.618, 0, -1, 3, 1.618, -0.618, 0, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Triply degenerate gerade representation.' },
      { symbol: 'Gg', id: 'Gg', dimension: 4, characters: [4, -1, -1, 1, 0, 4, -1, -1, 1, 0], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Four-fold degenerate gerade representation.' },
      { symbol: 'Hg', id: 'Hg', dimension: 5, characters: [5, 0, 0, -1, 1, 5, 0, 0, -1, 1], linearFunctions: [], quadraticFunctions: ['(2z²-x²-y², x²-y², xy, xz, yz)'], orbitalBasis: 'All five d-orbitals (C60 hu/hg)', irActive: false, ramanActive: true, description: 'Five-fold degenerate gerade representation spanning all five d-orbitals.' },
      { symbol: 'Au', id: 'Au', dimension: 1, characters: [1, 1, 1, 1, 1, -1, -1, -1, -1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: '', irActive: false, ramanActive: false, description: 'Singly degenerate ungerade representation.' },
      { symbol: 'T1u', id: 'T1u', dimension: 3, characters: [3, 1.618, -0.618, 0, -1, -3, 0.618, -1.618, 0, 1], linearFunctions: ['(x, y, z)'], quadraticFunctions: [], orbitalBasis: '(px, py, pz) (C60 t1u LUMO)', irActive: true, ramanActive: false, description: 'Triply degenerate IR-active representation (C60 LUMO).' },
      { symbol: 'T2u', id: 'T2u', dimension: 3, characters: [3, -0.618, 1.618, 0, -1, -3, -1.618, 0.618, 0, 1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'f-orbitals', irActive: false, ramanActive: false, description: 'Triply degenerate ungerade f-orbital representation.' },
      { symbol: 'Gu', id: 'Gu', dimension: 4, characters: [4, -1, -1, 1, 0, -4, 1, 1, -1, 0], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'f-orbitals', irActive: false, ramanActive: false, description: 'Four-fold degenerate ungerade representation.' },
      { symbol: 'Hu', id: 'Hu', dimension: 5, characters: [5, 0, 0, -1, 1, -5, 0, 0, 1, -1], linearFunctions: [], quadraticFunctions: [], orbitalBasis: 'C60 hu HOMO', irActive: false, ramanActive: false, description: 'Five-fold degenerate ungerade representation (HOMO of C60).' },
    ],
    elements3D: [
      { id: 'c5_1', type: 'axis', label: 'C5 (Icosahedral)', schoenflies: 'C5', fold: 5, axis: [0, 0.5257, 0.8507], color: '#2563EB', description: '5-fold icosahedral axis' },
      { id: 'c3_1', type: 'axis', label: 'C3 (Icosahedral)', schoenflies: 'C3', fold: 3, axis: [0.577, 0.577, 0.577], color: '#0284C7', description: '3-fold icosahedral axis' },
      { id: 'inv', type: 'center', label: 'Inversion Center i', schoenflies: 'i', color: '#DC2626', description: 'Center of inversion' },
    ],
    defaultReduciblePresets: [
      { id: 'ih_12', label: 'Γ_12 (12 Vertices of Icosahedron [B12H12]2-)', characters: [12, 2, 2, 0, 0, 0, 0, 0, 0, 4], description: 'Decomposes into Ag + T1u + Hg + T2u' },
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

    const lowerSymbol = irrep.symbol.toLowerCase().replace(/ς/g, 'σ');

    results.push({
      irrepSymbol: lowerSymbol,
      irrepId: irrep.id,
      coefficient,
      sumContribution,
      terms,
    });

    if (coefficient > 0.05) {
      const roundedInt = Math.abs(coefficient - Math.round(coefficient)) < 0.15
        ? Math.round(coefficient)
        : coefficient;
      formulaParts.push(roundedInt === 1 ? lowerSymbol : `${roundedInt}${lowerSymbol}`);
    }
  }

  return {
    formula: formulaParts.length > 0 ? formulaParts.join(' + ') : '0',
    results,
    isValid: formulaParts.length > 0,
  };
}
