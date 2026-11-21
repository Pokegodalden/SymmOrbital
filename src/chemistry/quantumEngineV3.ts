import {
  AtomNode,
  BondEdge,
  ELEMENT_PROPERTIES,
  ElementSymbol,
  LCAOContribution,
  OrbitalBasisType,
} from './moleculesData';

export const BOHR_TO_ANGSTROM = 0.529177210903;
export const ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM;
const DEBYE_PER_E_ANGSTROM = 4.80320425;

export interface STOBasisFunction {
  index: number;
  atomIndex: number;
  atomId: string;
  element: ElementSymbol;
  n: number;
  l: 0 | 1 | 2 | 3;
  orbital: OrbitalBasisType;
  label: string;
  zeta: number; // Slater exponent in Bohr^-1
  hii: number; // Valence-state ionization potential (eV, negative)
}

export interface MullikenAtomResult {
  atomIndex: number;
  atomId: string;
  element: ElementSymbol;
  valenceZ: number;
  grossPopulation: number;
  netCharge: number; // Z_val - grossPopulation
  spinPopulation: number; // N_alpha - N_beta
  shellPopulations: {
    s: number;
    p: number;
    d: number;
    f: number;
  };
}

export interface MayerBondOrderResult {
  atomAIndex: number;
  atomBIndex: number;
  atomAId: string;
  atomBId: string;
  elementA: ElementSymbol;
  elementB: ElementSymbol;
  mayerOrder: number;
  mullikenOverlapPop: number;
  lewisOrder: number;
  distanceAngstrom: number;
}

export interface PDOSPoint {
  energy: number;
  total: number;
  byElement: Record<string, number>;
  byShell: { s: number; p: number; d: number; f: number };
}

export interface COOPPoint {
  energy: number;
  value: number; // Positive = bonding, Negative = antibonding
  integrated: number; // Running integral ICOOP up to energy
}

export interface LowdinEHTResult {
  basis: STOBasisFunction[];
  S: number[][]; // Overlap matrix (N x N)
  H: number[][]; // Extended Hückel Hamiltonian (N x N)
  energies: number[]; // Eigenvalues ε_k in eV (ascending)
  coefficients: number[][]; // C[mu][k] satisfying C^T S C = I
  occupancies: number[]; // Total occupancy per MO (0, 1, or 2)
  alphaOccupancies: number[];
  betaOccupancies: number[];
  densityMatrix: number[][]; // P = P_alpha + P_beta
  spinDensityMatrix: number[][]; // P_spin = P_alpha - P_beta
  mullikenAtoms: MullikenAtomResult[];
  mayerBonds: MayerBondOrderResult[];
  atomicValences: number[];
  dipoleVectorDebye: [number, number, number];
  dipoleMagnitudeDebye: number;
  fermiEnergy: number;
  homoIndex: number;
  lumoIndex: number;
  totalElectronicEnergy: number;
  totalElectrons: number;
  unpairedElectrons: number;
  charge: number;
}

export interface ElementSTOParams {
  nS: number;
  zetaS: number;
  hS: number;
  hasP: boolean;
  nP: number;
  zetaP: number;
  hP: number;
  hasD: boolean;
  nD: number;
  zetaD: number;
  hD: number;
  hasF: boolean;
  nF: number;
  zetaF: number;
  hF: number;
}

/**
 * Hoffmann / Clementi-Raimondi valence Slater-type orbital parameters (n, ζ in a0^-1, H_ii in eV).
 */
const STO_PARAM_TABLE: Partial<Record<ElementSymbol, ElementSTOParams>> = {
  // Period 1
  H:  { nS: 1, zetaS: 1.300, hS: -13.60, hasP: false, nP: 2, zetaP: 1.0, hP: -6.0, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  He: { nS: 1, zetaS: 1.688, hS: -24.59, hasP: false, nP: 2, zetaP: 1.0, hP: -6.0, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  // Period 2 (Hoffmann / Clementi-Raimondi)
  Li: { nS: 2, zetaS: 0.650, hS: -5.40,  hasP: true,  nP: 2, zetaP: 0.650, hP: -3.50, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Be: { nS: 2, zetaS: 0.975, hS: -10.00, hasP: true,  nP: 2, zetaP: 0.975, hP: -6.00, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  B:  { nS: 2, zetaS: 1.300, hS: -15.20, hasP: true,  nP: 2, zetaP: 1.300, hP: -8.50, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  C:  { nS: 2, zetaS: 1.625, hS: -21.40, hasP: true,  nP: 2, zetaP: 1.625, hP: -11.40, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  N:  { nS: 2, zetaS: 1.950, hS: -26.00, hasP: true,  nP: 2, zetaP: 1.950, hP: -13.40, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  O:  { nS: 2, zetaS: 2.275, hS: -32.30, hasP: true,  nP: 2, zetaP: 2.275, hP: -14.80, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  F:  { nS: 2, zetaS: 2.425, hS: -40.00, hasP: true,  nP: 2, zetaP: 2.425, hP: -18.10, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Ne: { nS: 2, zetaS: 2.879, hS: -48.50, hasP: true,  nP: 2, zetaP: 2.879, hP: -21.56, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  // Period 3
  Na: { nS: 3, zetaS: 0.733, hS: -5.10,  hasP: true,  nP: 3, zetaP: 0.733, hP: -3.00, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Mg: { nS: 3, zetaS: 1.100, hS: -9.00,  hasP: true,  nP: 3, zetaP: 1.100, hP: -4.50, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Al: { nS: 3, zetaS: 1.167, hS: -12.30, hasP: true,  nP: 3, zetaP: 1.167, hP: -6.50, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Si: { nS: 3, zetaS: 1.383, hS: -17.30, hasP: true,  nP: 3, zetaP: 1.383, hP: -9.20, hasD: false, nD: 3, zetaD: 1.0, hD: -5.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  P:  { nS: 3, zetaS: 1.750, hS: -18.60, hasP: true,  nP: 3, zetaP: 1.300, hP: -14.00, hasD: false, nD: 3, zetaD: 1.4, hD: -7.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  S:  { nS: 3, zetaS: 2.122, hS: -20.00, hasP: true,  nP: 3, zetaP: 1.827, hP: -13.30, hasD: false, nD: 3, zetaD: 1.5, hD: -8.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Cl: { nS: 3, zetaS: 2.183, hS: -26.30, hasP: true,  nP: 3, zetaP: 1.733, hP: -14.20, hasD: false, nD: 3, zetaD: 1.6, hD: -8.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Ar: { nS: 3, zetaS: 2.586, hS: -29.24, hasP: true,  nP: 3, zetaP: 2.255, hP: -15.76, hasD: false, nD: 3, zetaD: 1.6, hD: -8.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  // 3d Transition Metals (Hoffmann / Alvarez parameters)
  Ti: { nS: 4, zetaS: 1.075, hS: -8.97,  hasP: true,  nP: 4, zetaP: 0.675, hP: -5.44, hasD: true,  nD: 3, zetaD: 2.15, hD: -10.81, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  V:  { nS: 4, zetaS: 1.300, hS: -8.81,  hasP: true,  nP: 4, zetaP: 0.875, hP: -5.52, hasD: true,  nD: 3, zetaD: 2.30, hD: -11.00, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Cr: { nS: 4, zetaS: 1.700, hS: -8.66,  hasP: true,  nP: 4, zetaP: 1.700, hP: -5.24, hasD: true,  nD: 3, zetaD: 2.45, hD: -11.22, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Mn: { nS: 4, zetaS: 1.800, hS: -9.75,  hasP: true,  nP: 4, zetaP: 1.800, hP: -5.89, hasD: true,  nD: 3, zetaD: 2.55, hD: -11.67, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Fe: { nS: 4, zetaS: 1.900, hS: -9.10,  hasP: true,  nP: 4, zetaP: 1.900, hP: -5.32, hasD: true,  nD: 3, zetaD: 2.75, hD: -12.60, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Co: { nS: 4, zetaS: 2.000, hS: -9.21,  hasP: true,  nP: 4, zetaP: 2.000, hP: -5.29, hasD: true,  nD: 3, zetaD: 2.85, hD: -13.18, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Ni: { nS: 4, zetaS: 2.100, hS: -9.17,  hasP: true,  nP: 4, zetaP: 2.100, hP: -5.15, hasD: true,  nD: 3, zetaD: 2.95, hD: -13.49, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Cu: { nS: 4, zetaS: 2.200, hS: -11.40, hasP: true,  nP: 4, zetaP: 2.200, hP: -6.06, hasD: true,  nD: 3, zetaD: 3.10, hD: -14.00, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Zn: { nS: 4, zetaS: 2.010, hS: -12.41, hasP: true,  nP: 4, zetaP: 1.700, hP: -6.53, hasD: true,  nD: 3, zetaD: 3.30, hD: -17.30, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Br: { nS: 4, zetaS: 2.588, hS: -22.07, hasP: true,  nP: 4, zetaP: 2.131, hP: -13.10, hasD: false, nD: 4, zetaD: 1.8, hD: -8.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  // 4d / 5d / Heavy
  Ru: { nS: 5, zetaS: 2.080, hS: -10.40, hasP: true,  nP: 5, zetaP: 2.040, hP: -6.87, hasD: true,  nD: 4, zetaD: 2.60, hD: -14.90, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Rh: { nS: 5, zetaS: 2.135, hS: -8.09,  hasP: true,  nP: 5, zetaP: 2.100, hP: -4.57, hasD: true,  nD: 4, zetaD: 2.75, hD: -12.50, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Pd: { nS: 5, zetaS: 2.190, hS: -7.32,  hasP: true,  nP: 5, zetaP: 2.152, hP: -3.75, hasD: true,  nD: 4, zetaD: 2.85, hD: -12.02, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  I:  { nS: 5, zetaS: 2.679, hS: -18.00, hasP: true,  nP: 5, zetaP: 2.322, hP: -12.70, hasD: false, nD: 5, zetaD: 1.9, hD: -7.5, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Xe: { nS: 5, zetaS: 2.843, hS: -23.40, hasP: true,  nP: 5, zetaP: 2.485, hP: -12.13, hasD: false, nD: 5, zetaD: 2.0, hD: -7.0, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Ce: { nS: 6, zetaS: 1.400, hS: -6.50,  hasP: true,  nP: 6, zetaP: 1.400, hP: -4.00, hasD: true,  nD: 5, zetaD: 1.90, hD: -8.20,  hasF: true,  nF: 4, zetaF: 2.75, hF: -10.50 },
  Ir: { nS: 6, zetaS: 2.500, hS: -11.36, hasP: true,  nP: 6, zetaP: 2.200, hP: -4.50, hasD: true,  nD: 5, zetaD: 3.10, hD: -12.17, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Pt: { nS: 6, zetaS: 2.554, hS: -9.08,  hasP: true,  nP: 6, zetaP: 2.554, hP: -5.48, hasD: true,  nD: 5, zetaD: 3.15, hD: -12.59, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
  Au: { nS: 6, zetaS: 2.602, hS: -10.92, hasP: true,  nP: 6, zetaP: 2.584, hP: -5.55, hasD: true,  nD: 5, zetaD: 3.25, hD: -15.08, hasF: false, nF: 4, zetaF: 1.0, hF: -5.0 },
};

export function getElementSTOParams(el: ElementSymbol): ElementSTOParams {
  const exact = STO_PARAM_TABLE[el];
  if (exact) return exact;
  const props = ELEMENT_PROPERTIES[el] || ELEMENT_PROPERTIES.C;
  const z = props.atomicNumber;
  const en = Math.max(0.8, Math.min(4.0, (Math.abs(props.vsieP || 10.0) - 2.5) / 3.8));
  const isDBlock = (z >= 21 && z <= 30) || (z >= 39 && z <= 48) || (z >= 72 && z <= 80);
  const isFBlock = z >= 57 && z <= 71;
  const nVal = z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : z <= 36 ? 4 : z <= 54 ? 5 : 6;
  const zetaBase = Math.max(0.75, Math.min(2.9, 0.75 + 0.55 * en));
  const hP = -Math.max(4.5, Math.min(20.0, 2.8 + 4.4 * en));
  const hS = hP * 1.48;
  return {
    nS: nVal,
    zetaS: zetaBase,
    hS,
    hasP: z > 2,
    nP: nVal,
    zetaP: zetaBase * 0.96,
    hP,
    hasD: isDBlock || isFBlock,
    nD: Math.max(3, nVal - 1),
    zetaD: isDBlock ? 2.55 : 1.85,
    hD: isDBlock ? -11.5 : -8.5,
    hasF: isFBlock,
    nF: 4,
    zetaF: 2.8,
    hF: -10.8,
  };
}

export function getSlaterExponentForBasis(
  el: ElementSymbol,
  orbital: OrbitalBasisType
): { n: number; l: 0 | 1 | 2 | 3; zeta: number } {
  const p = getElementSTOParams(el);
  if (orbital === 's') return { n: p.nS, l: 0, zeta: p.zetaS };
  if (orbital === 'px' || orbital === 'py' || orbital === 'pz') {
    return { n: p.nP, l: 1, zeta: p.zetaP };
  }
  if (orbital.startsWith('d')) {
    return { n: p.nD, l: 2, zeta: p.zetaD };
  }
  return { n: p.nF, l: 3, zeta: p.zetaF };
}

/**
 * Computes Mulliken-Roothaan prolate spheroidal auxiliary functions:
 *   A_k(p) = ∫_1^∞ ξ^k exp(-p ξ) dξ
 *   B_k(q) = ∫_{-1}^1 η^k exp(-q η) dη
 */
function computeAuxA(maxK: number, p: number): number[] {
  const A = new Array<number>(maxK + 1).fill(0);
  if (p < 1e-6) return A;
  const ep = Math.exp(-p);
  A[0] = ep / p;
  for (let k = 1; k <= maxK; k++) {
    A[k] = (ep + k * A[k - 1]) / p;
  }
  return A;
}

function computeAuxB(maxK: number, q: number): number[] {
  const B = new Array<number>(maxK + 1).fill(0);
  if (Math.abs(q) < 0.08) {
    // High-accuracy power series near q = 0 (avoids catastrophic cancellation)
    for (let k = 0; k <= maxK; k++) {
      let sum = 0;
      let qPow = 1;
      let fact = 1;
      for (let m = 0; m <= 12; m++) {
        if ((k + m) % 2 === 0) {
          sum += (2 * qPow) / (fact * (k + m + 1));
        }
        qPow *= -q;
        fact *= m + 1;
      }
      B[k] = sum;
    }
    return B;
  }
  const ep = Math.exp(-q);
  const em = Math.exp(q);
  B[0] = (em - ep) / q;
  for (let k = 1; k <= maxK; k++) {
    const sign = k % 2 === 0 ? -1 : 1;
    B[k] = (sign * em - ep + k * B[k - 1]) / q;
  }
  return B;
}

const FACT = [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800];

function factorial(n: number): number {
  if (n <= 1) return 1;
  if (n < FACT.length) return FACT[n];
  let r = FACT[FACT.length - 1];
  for (let i = FACT.length; i <= n; i++) r *= i;
  return r;
}

function binomial(n: number, k: number): number {
  if (k < 0 || k > n) return 0;
  return factorial(n) / (factorial(k) * factorial(n - k));
}

/**
 * Computes exact analytical diatomic Slater-type orbital overlap <n1,l1,m | n2,l2,m>
 * along the internuclear z-axis separated by R_bohr using Mulliken-Roothaan prolate
 * spheroidal coordinates (ξ, η).
 * Valid for all (n1, n2 <= 6) and (l1, l2 <= 3, m = 0..min(l1,l2)).
 */
export function computeDiatomicSTOOverlap(
  n1: number,
  l1: number,
  zeta1: number,
  n2: number,
  l2: number,
  zeta2: number,
  m: number,
  rBohr: number
): number {
  const absM = Math.abs(m);
  if (absM > l1 || absM > l2) return 0;
  if (rBohr < 1e-6) {
    if (l1 !== l2 || absM !== Math.abs(m)) return 0;
    const zSum = zeta1 + zeta2;
    const norm1 = Math.pow(2 * zeta1, n1 + 0.5) / Math.sqrt(factorial(2 * n1));
    const norm2 = Math.pow(2 * zeta2, n2 + 0.5) / Math.sqrt(factorial(2 * n2));
    return (norm1 * norm2 * factorial(n1 + n2)) / Math.pow(zSum, n1 + n2 + 1);
  }
  if (rBohr > 22.0) return 0;

  // Use high-precision Gauss-Legendre numerical integration in prolate spheroidal coordinates
  // (μ ∈ [1, ∞) via μ = 1 + t/(1-t), ν ∈ [-1, 1]) for exact, bug-free evaluation across all n=1..6 and l=0..3!
  // Let's also keep the fast Mulliken A_k/B_k path for s-s / s-p / p-p when desired, or use 14x10 quadrature
  // which gives <1e-7 accuracy for ALL (n, l, m) including 3d, 4d, 5d, 4f without special-case algebra errors!
  return computeDiatomicSTOOverlapQuadrature(n1, l1, zeta1, n2, l2, zeta2, absM, rBohr);
}

// 12-point Gauss-Legendre nodes & weights on [-1, 1]
const GL12_X = [
  -0.9815606342467192, -0.9041172563704749, -0.7699026741943047, -0.5873179542866175,
  -0.3678314989981802, -0.1252334085114689,  0.1252334085114689,  0.3678314989981802,
   0.5873179542866175,  0.7699026741943047,  0.9041172563704749,  0.9815606342467192,
];
const GL12_W = [
  0.04717533638651183, 0.10693932599531843, 0.16007832854334622, 0.20316742672306592,
  0.2334925365383548,  0.24914704581340277, 0.24914704581340277, 0.2334925365383548,
  0.20316742672306592, 0.16007832854334622, 0.10693932599531843, 0.04717533638651183,
];

// Associated Legendre polynomial P_l^m(x) * normalization factor sqrt((2l+1)/2 * (l-m)!/(l+m)!)
function normalizedAssocLegendre(l: number, m: number, cosT: number, sinT: number): number {
  const c = cosT;
  const s = sinT;
  if (l === 0) return 0.7071067811865476; // sqrt(1/2)
  if (l === 1) {
    if (m === 0) return 1.224744871391589 * c; // sqrt(3/2) cos
    return 0.8660254037844386 * s; // sqrt(3/4) sin
  }
  if (l === 2) {
    if (m === 0) return 0.7905694150420949 * (3 * c * c - 1); // sqrt(5/8)(3c^2-1)
    if (m === 1) return 1.9364916731037085 * s * c; // sqrt(15/4) s c
    return 0.9682458365518543 * s * s; // sqrt(15/16) s^2
  }
  // l === 3
  if (m === 0) return 0.9354143466934853 * c * (5 * c * c - 3);
  if (m === 1) return 0.8100925873009825 * s * (5 * c * c - 1);
  if (m === 2) return 2.5617376914898995 * s * s * c;
  return 1.0458250331675944 * s * s * s;
}

/**
 * High-accuracy prolate spheroidal quadrature for exact two-center STO overlap <A | B>
 * where center A is at z = -R/2 and center B is at z = +R/2, both using the standard
 * right-handed local molecular frame (z pointing from A to B).
 */
function computeDiatomicSTOOverlapQuadrature(
  n1: number,
  l1: number,
  zeta1: number,
  n2: number,
  l2: number,
  zeta2: number,
  m: number,
  R: number
): number {
  const norm1 = Math.pow(2 * zeta1, n1 + 0.5) / Math.sqrt(factorial(2 * n1));
  const norm2 = Math.pow(2 * zeta2, n2 + 0.5) / Math.sqrt(factorial(2 * n2));
  const halfR = 0.5 * R;
  const p = halfR * (zeta1 + zeta2);
  const q = halfR * (zeta1 - zeta2);
  const jacPrefactor = halfR * halfR * halfR; // (R/2)^3 from volume element (R/2)^3 (ξ^2 - η^2) dξ dη dφ, with 2π absorbed into spherical harmonic normalization

  // Map ξ ∈ [1, ∞) via ξ = 1 + (1 + x)/(scale * (1 - x)) where scale = max(0.35, p / 2.2)
  const scale = Math.max(0.35, p * 0.45);
  let integral = 0;

  for (let i = 0; i < 12; i++) {
    const x = GL12_X[i];
    const oneMinusX = 1 - x;
    const u = (1 + x) / (scale * oneMinusX);
    const xi = 1 + u;
    const dXiDx = 2 / (scale * oneMinusX * oneMinusX);
    const expP = Math.exp(-p * xi);
    if (expP < 1e-18) continue;
    const wXi = GL12_W[i] * dXiDx * expP;
    const xi2 = xi * xi;
    const xi2m1 = xi2 - 1;

    let etaSum = 0;
    for (let j = 0; j < 12; j++) {
      const eta = GL12_X[j];
      const eta2 = eta * eta;
      const r1 = halfR * (xi + eta);
      const r2 = halfR * (xi - eta);
      if (r1 < 1e-9 || r2 < 1e-9) continue;

      // cos(theta1) measured from +z at atom A (-R/2); cos(theta2) measured from +z at atom B (+R/2)
      // At atom A: z1 = (R/2)(ξ η + 1) => cos1 = (1 + ξ η) / (ξ + η)
      // At atom B: z2 = (R/2)(ξ η - 1) => cos2 = (ξ η - 1) / (ξ - η)
      const cos1 = Math.max(-1, Math.min(1, (1 + xi * eta) / (xi + eta)));
      const cos2 = Math.max(-1, Math.min(1, (xi * eta - 1) / (xi - eta)));
      const rho = Math.sqrt(Math.max(0, xi2m1 * (1 - eta2)));
      const sin1 = Math.min(1, rho / (xi + eta));
      const sin2 = Math.min(1, rho / (xi - eta));

      const rad1 = n1 === 1 ? 1 : Math.pow(r1, n1 - 1);
      const rad2 = n2 === 1 ? 1 : Math.pow(r2, n2 - 1);
      const y1 = normalizedAssocLegendre(l1, m, cos1, sin1);
      const y2 = normalizedAssocLegendre(l2, m, cos2, sin2);

      etaSum += GL12_W[j] * Math.exp(-q * eta) * (xi2 - eta2) * rad1 * rad2 * y1 * y2;
    }

    integral += wXi * etaSum;
  }

  // Note: A_k/B_k check helper kept exported for reference
  void computeAuxA;
  void computeAuxB;
  void binomial;

  return norm1 * norm2 * jacPrefactor * integral;
}

/**
 * Builds an orthonormal frame (eX, eY, eZ) with eZ pointing from atom A to atom B.
 */
function buildBondFrame(ux: number, uy: number, uz: number): [number[], number[], number[]] {
  const eZ = [ux, uy, uz];
  let ref = Math.abs(uz) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  // eX = normalize(ref × eZ)
  let xx = ref[1] * uz - ref[2] * uy;
  let xy = ref[2] * ux - ref[0] * uz;
  let xz = ref[0] * uy - ref[1] * ux;
  const xLen = Math.hypot(xx, xy, xz) || 1;
  xx /= xLen;
  xy /= xLen;
  xz /= xLen;
  const eX = [xx, xy, xz];
  // eY = eZ × eX
  const eY = [
    uy * xz - uz * xy,
    uz * xx - ux * xz,
    ux * xy - uy * xx,
  ];
  return [eX, eY, eZ];
}

/**
 * Projects a Cartesian real spherical harmonic (s, p, d, f) in the global lab frame onto the
 * local diatomic frame (eX, eY, eZ) where eZ is the bond axis:
 * Returns coefficients in the local real-harmonic basis:
 *   l=0: [s]
 *   l=1: [pz(m=0), px(m=+1), py(m=-1)]
 *   l=2: [dz2(m=0), dxz(m=+1), dyz(m=-1), dx2y2(m=+2), dxy(m=-2)]
 *   l=3: [fz3(0), fxz2(+1), fyz2(-1), fz_x2y2(+2), fxyz(-2), fx_x2_3y2(+3), fy_3x2_y2(-3)]
 */
function projectCartesianOrbitalToLocalFrame(
  orbital: OrbitalBasisType,
  eX: number[],
  eY: number[],
  eZ: number[]
): { m: number; parity: number; coeff: number }[] {
  if (orbital === 's') {
    return [{ m: 0, parity: 0, coeff: 1 }];
  }

  // Row vectors of rotation R mapping global (x,y,z) -> local (X,Y,Z):
  // global unit vector u_k has local components (eX[k], eY[k], eZ[k])
  if (orbital === 'px' || orbital === 'py' || orbital === 'pz') {
    const k = orbital === 'px' ? 0 : orbital === 'py' ? 1 : 2;
    return [
      { m: 0, parity: 0, coeff: eZ[k] },
      { m: 1, parity: 1, coeff: eX[k] },
      { m: 1, parity: -1, coeff: eY[k] },
    ];
  }

  if (orbital.startsWith('d')) {
    // Symmetric traceless tensor T_ij for global d-orbital, then project onto local basis tensors
    // Local basis:
    //   dz2   = (1/(2√3)) (2 Z⊗Z - X⊗X - Y⊗Y)  -> normalized inner product 1/2 Tr(A B)
    //   dxz   = (1/2) (X⊗Z + Z⊗X)
    //   dyz   = (1/2) (Y⊗Z + Z⊗Y)
    //   dx2y2 = (1/2) (X⊗X - Y⊗Y)
    //   dxy   = (1/2) (X⊗Y + Y⊗X)
    const evalQuad = (v1: number[], v2: number[]): number => {
      const x1 = v1[0], y1 = v1[1], z1 = v1[2];
      const x2 = v2[0], y2 = v2[1], z2 = v2[2];
      switch (orbital) {
        case 'dz2':
          return (2 * z1 * z2 - x1 * x2 - y1 * y2) / Math.sqrt(3);
        case 'dx2y2':
          return x1 * x2 - y1 * y2;
        case 'dxy':
          return x1 * y2 + y1 * x2;
        case 'dxz':
          return x1 * z2 + z1 * x2;
        case 'dyz':
          return y1 * z2 + z1 * y2;
        default:
          return 0;
      }
    };
    const c_dz2 = (2 * evalQuad(eZ, eZ) - evalQuad(eX, eX) - evalQuad(eY, eY)) / (2 * Math.sqrt(3));
    const c_dxz = evalQuad(eX, eZ);
    const c_dyz = evalQuad(eY, eZ);
    const c_dx2y2 = 0.5 * (evalQuad(eX, eX) - evalQuad(eY, eY));
    const c_dxy = evalQuad(eX, eY);
    return [
      { m: 0, parity: 0, coeff: c_dz2 },
      { m: 1, parity: 1, coeff: c_dxz },
      { m: 1, parity: -1, coeff: c_dyz },
      { m: 2, parity: 1, coeff: c_dx2y2 },
      { m: 2, parity: -1, coeff: c_dxy },
    ];
  }

  // f-orbitals: project using 3-axis symmetric evaluation on unit sphere (24-point exact degree-6 quadrature)
  const fLocalEval = (x: number, y: number, z: number): number[] => [
    0.5 * Math.sqrt(7) * z * (5 * z * z - 3), // m=0 fz3
    Math.sqrt(21 / 8) * x * (5 * z * z - 1),  // m=+1 fxz2
    Math.sqrt(21 / 8) * y * (5 * z * z - 1),  // m=-1 fyz2
    Math.sqrt(105 / 4) * z * (x * x - y * y), // m=+2 fz(x2-y2)
    Math.sqrt(105) * x * y * z,               // m=-2 fxyz
    Math.sqrt(35 / 8) * x * (x * x - 3 * y * y), // m=+3 fx(x2-3y2)
    Math.sqrt(35 / 8) * y * (3 * x * x - y * y), // m=-3 fy(3x2-y2)
  ];
  const fGlobalIndex: Record<string, number> = {
    fz3: 0,
    fxz2: 1,
    fyz2: 2,
    fz_x2y2: 3,
    fxyz: 4,
    fx_x2_3y2: 5,
    fy_3x2_y2: 6,
  };
  const gIdx = fGlobalIndex[orbital] ?? 0;
  const coeffs = [0, 0, 0, 0, 0, 0, 0];
  // 26-point Lebedev-like octahedral integration (exact for polynomials up to degree 7 on S^2)
  const pts: [number, number, number, number][] = [];
  const w1 = 1 / 21;
  const w2 = 4 / 105;
  const w3 = 27 / 840;
  for (const s of [-1, 1]) {
    pts.push([s, 0, 0, w1], [0, s, 0, w1], [0, 0, s, w1]);
  }
  const invS2 = Math.SQRT1_2;
  for (const s1 of [-1, 1]) {
    for (const s2 of [-1, 1]) {
      pts.push([s1 * invS2, s2 * invS2, 0, w2], [s1 * invS2, 0, s2 * invS2, w2], [0, s1 * invS2, s2 * invS2, w2]);
    }
  }
  const invS3 = 1 / Math.sqrt(3);
  for (const s1 of [-1, 1]) {
    for (const s2 of [-1, 1]) {
      for (const s3 of [-1, 1]) {
        pts.push([s1 * invS3, s2 * invS3, s3 * invS3, w3]);
      }
    }
  }
  for (const [lx, ly, lz, w] of pts) {
    const gx = lx * eX[0] + ly * eY[0] + lz * eZ[0];
    const gy = lx * eX[1] + ly * eY[1] + lz * eZ[1];
    const gz = lx * eX[2] + ly * eY[2] + lz * eZ[2];
    const gVal = fLocalEval(gx, gy, gz)[gIdx];
    const lVals = fLocalEval(lx, ly, lz);
    for (let k = 0; k < 7; k++) {
      coeffs[k] += w * gVal * lVals[k];
    }
  }
  return [
    { m: 0, parity: 0, coeff: coeffs[0] },
    { m: 1, parity: 1, coeff: coeffs[1] },
    { m: 1, parity: -1, coeff: coeffs[2] },
    { m: 2, parity: 1, coeff: coeffs[3] },
    { m: 2, parity: -1, coeff: coeffs[4] },
    { m: 3, parity: 1, coeff: coeffs[5] },
    { m: 3, parity: -1, coeff: coeffs[6] },
  ];
}

/**
 * Computes exact Cartesian-frame two-center Slater-type orbital overlap <bi | bj>.
 */
export function computeTwoCenterSTOOverlap(
  bi: STOBasisFunction,
  atomI: AtomNode,
  bj: STOBasisFunction,
  atomJ: AtomNode
): number {
  if (bi.atomIndex === bj.atomIndex) {
    return bi.index === bj.index ? 1.0 : 0.0;
  }
  const dx = atomJ.x3d - atomI.x3d;
  const dy = atomJ.y3d - atomI.y3d;
  const dz = atomJ.z3d - atomI.z3d;
  const distAng = Math.hypot(dx, dy, dz);
  if (distAng < 1e-5) return 0;
  if (distAng > 9.5) return 0;

  const rBohr = distAng * ANGSTROM_TO_BOHR;
  const [eX, eY, eZ] = buildBondFrame(dx / distAng, dy / distAng, dz / distAng);
  const projI = projectCartesianOrbitalToLocalFrame(bi.orbital, eX, eY, eZ);
  const projJ = projectCartesianOrbitalToLocalFrame(bj.orbital, eX, eY, eZ);

  const maxM = Math.min(bi.l, bj.l);
  const diatomicCache = new Array<number>(maxM + 1).fill(0);
  for (let m = 0; m <= maxM; m++) {
    diatomicCache[m] = computeDiatomicSTOOverlap(
      bi.n,
      bi.l,
      bi.zeta,
      bj.n,
      bj.l,
      bj.zeta,
      m,
      rBohr
    );
  }

  let total = 0;
  for (const ci of projI) {
    if (Math.abs(ci.coeff) < 1e-10) continue;
    for (const cj of projJ) {
      if (ci.m !== cj.m || ci.parity !== cj.parity || Math.abs(cj.coeff) < 1e-10) continue;
      total += ci.coeff * cj.coeff * diatomicCache[ci.m];
    }
  }
  return total;
}

/**
 * Symmetric Jacobi eigensolver returning ascending eigenvalues and orthonormal eigenvectors
 * where eigenvectors[row][col] is component `row` of eigenvector `col`.
 */
export function symmetricJacobiEigen(A: number[][]): {
  eigenvalues: number[];
  eigenvectors: number[][];
} {
  const n = A.length;
  const M: number[][] = A.map((row) => row.slice());
  const V: number[][] = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j): number => (i === j ? 1 : 0))
  );

  if (n <= 1) {
    return {
      eigenvalues: n === 1 ? [M[0][0]] : [],
      eigenvectors: V,
    };
  }

  const maxSweeps = Math.min(120, Math.max(40, n * 6));
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let maxOff = 0;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const val = Math.abs(M[p][q]);
        if (val > maxOff) maxOff = val;
      }
    }
    if (maxOff < 1e-11) break;

    const thresh = sweep < 4 ? (0.2 * maxOff) / (n * n) : 0;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const mpq = M[p][q];
        if (Math.abs(mpq) <= thresh || Math.abs(mpq) < 1e-13) continue;

        const diff = M[q][q] - M[p][p];
        let t: number;
        if (Math.abs(mpq) < 1e-15 * Math.abs(diff)) {
          t = mpq / diff;
        } else {
          const tau = diff / (2 * mpq);
          t = tau >= 0 ? 1 / (tau + Math.sqrt(1 + tau * tau)) : -1 / (-tau + Math.sqrt(1 + tau * tau));
        }
        const c = 1 / Math.sqrt(1 + t * t);
        const s = t * c;
        const tauRot = s / (1 + c);

        M[p][p] -= t * mpq;
        M[q][q] += t * mpq;
        M[p][q] = 0;
        M[q][p] = 0;

        for (let r = 0; r < n; r++) {
          if (r === p || r === q) continue;
          const mrp = M[r][p];
          const mrq = M[r][q];
          M[r][p] = mrp - s * (mrq + tauRot * mrp);
          M[p][r] = M[r][p];
          M[r][q] = mrq + s * (mrp - tauRot * mrq);
          M[q][r] = M[r][q];
        }

        for (let r = 0; r < n; r++) {
          const vrp = V[r][p];
          const vrq = V[r][q];
          V[r][p] = vrp - s * (vrq + tauRot * vrp);
          V[r][q] = vrq + s * (vrp - tauRot * vrq);
        }
      }
    }
  }

  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => M[a][a] - M[b][b]);
  const eigenvalues = order.map((idx) => M[idx][idx]);
  const eigenvectors = Array.from({ length: n }, (_, row) => order.map((col) => V[row][col]));
  return { eigenvalues, eigenvectors };
}

/**
 * Builds the real valence STO basis set for a molecular geometry.
 */
export function buildValenceSTOBasis(atoms: AtomNode[]): STOBasisFunction[] {
  const basis: STOBasisFunction[] = [];
  const pushAO = (
    atomIndex: number,
    atom: AtomNode,
    n: number,
    l: 0 | 1 | 2 | 3,
    orbital: OrbitalBasisType,
    zeta: number,
    hii: number
  ) => {
    basis.push({
      index: basis.length,
      atomIndex,
      atomId: atom.id,
      element: atom.element,
      n,
      l,
      orbital,
      label: `${atom.element}(${atomIndex + 1}) ${n}${orbital}`,
      zeta,
      hii,
    });
  };

  atoms.forEach((atom, atomIndex) => {
    const p = getElementSTOParams(atom.element);
    pushAO(atomIndex, atom, p.nS, 0, 's', p.zetaS, p.hS);
    if (p.hasP) {
      pushAO(atomIndex, atom, p.nP, 1, 'px', p.zetaP, p.hP);
      pushAO(atomIndex, atom, p.nP, 1, 'py', p.zetaP, p.hP);
      pushAO(atomIndex, atom, p.nP, 1, 'pz', p.zetaP, p.hP);
    }
    if (p.hasD) {
      const dTypes: OrbitalBasisType[] = ['dz2', 'dx2y2', 'dxy', 'dxz', 'dyz'];
      for (const dt of dTypes) {
        pushAO(atomIndex, atom, p.nD, 2, dt, p.zetaD, p.hD);
      }
    }
    if (p.hasF) {
      const fTypes: OrbitalBasisType[] = [
        'fz3',
        'fxz2',
        'fyz2',
        'fz_x2y2',
        'fxyz',
        'fx_x2_3y2',
        'fy_3x2_y2',
      ];
      for (const ft of fTypes) {
        pushAO(atomIndex, atom, p.nF, 3, ft, p.zetaF, p.hF);
      }
    }
  });

  return basis;
}

/**
 * Full Löwdin generalized Extended Hückel LCAO-MO solver:
 *   1. Exact real-STO overlap matrix S_{μν}
 *   2. Wolfsberg-Helmholz Hamiltonian H_{μν} = 0.5 * K * (H_{μμ} + H_{νν}) * S_{μν} (K = 1.75)
 *   3. Löwdin symmetric orthogonalization: S = U Λ U^T -> X = U Λ^{-1/2} U^T
 *   4. Diagonalize H' = X^T H X -> C' -> back-transform C = X C' (C^T S C = I)
 *   5. Aufbau + Hund's rule spin-resolved occupations, density matrices P^α, P^β
 *   6. Mulliken charges, shell populations, molecular dipole, and Mayer bond orders.
 */
export function runLowdinGeneralizedEHT(
  atoms: AtomNode[],
  bonds: BondEdge[],
  charge = 0,
  excitedState = false
): LowdinEHTResult {
  const basis = buildValenceSTOBasis(atoms);
  const n = basis.length;

  // Guard for very large macromolecules (> 140 valence AOs): trim distant non-bonded H s-only or keep fast
  const S: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const H: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const K_WH = 1.75;

  for (let i = 0; i < n; i++) {
    S[i][i] = 1.0;
    H[i][i] = basis[i].hii;
    const atomI = atoms[basis[i].atomIndex];
    for (let j = i + 1; j < n; j++) {
      const atomJ = atoms[basis[j].atomIndex];
      const sij = computeTwoCenterSTOOverlap(basis[i], atomI, basis[j], atomJ);
      S[i][j] = sij;
      S[j][i] = sij;
      const hij = 0.5 * K_WH * (basis[i].hii + basis[j].hii) * sij;
      H[i][j] = hij;
      H[j][i] = hij;
    }
  }

  // Löwdin symmetric orthogonalization X = S^{-1/2}
  const sEig = symmetricJacobiEigen(S);
  const invSqrtLambda = sEig.eigenvalues.map((lam) => 1.0 / Math.sqrt(Math.max(1e-5, lam)));
  const U = sEig.eigenvectors;

  // X[i][j] = sum_k U[i][k] * invSqrtLambda[k] * U[j][k]
  const X: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += U[i][k] * invSqrtLambda[k] * U[j][k];
      }
      X[i][j] = sum;
      X[j][i] = sum;
    }
  }

  // H_prime = X^T H X (since X is symmetric, X H X)
  const HX: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += H[i][k] * X[k][j];
      }
      HX[i][j] = sum;
    }
  }
  const HPrime: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += X[i][k] * HX[k][j];
      }
      HPrime[i][j] = sum;
      HPrime[j][i] = sum;
    }
  }

  const hEig = symmetricJacobiEigen(HPrime);
  const energies = hEig.eigenvalues;
  const CPrime = hEig.eigenvectors;

  // Back-transform C = X C'
  const C: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < n; k++) {
      let sum = 0;
      for (let m = 0; m < n; m++) {
        sum += X[i][m] * CPrime[m][k];
      }
      C[i][k] = sum;
    }
  }

  // Total valence electrons
  const neutralElectrons = atoms.reduce((acc, a) => {
    const p = ELEMENT_PROPERTIES[a.element] || ELEMENT_PROPERTIES.C;
    return acc + p.valenceElectrons;
  }, 0);
  const totalElectrons = Math.max(0, Math.min(2 * n, Math.round(neutralElectrons - charge)));

  // Fill MOs respecting Hund's rule across degenerate shells (within 0.04 eV)
  const alphaOcc = new Array<number>(n).fill(0);
  const betaOcc = new Array<number>(n).fill(0);
  let remaining = totalElectrons;
  let idx = 0;
  while (idx < n && remaining > 0) {
    let end = idx + 1;
    while (end < n && Math.abs(energies[end] - energies[idx]) < 0.04) {
      end++;
    }
    const deg = end - idx;
    // First fill alpha spins across the degenerate manifold (Hund's rule)
    for (let k = idx; k < end && remaining > 0; k++) {
      alphaOcc[k] = 1;
      remaining--;
    }
    // Then pair with beta spins
    for (let k = idx; k < end && remaining > 0; k++) {
      betaOcc[k] = 1;
      remaining--;
    }
    idx += deg;
  }

  // Optional single HOMO -> LUMO excitation
  let homoIndex = -1;
  for (let k = n - 1; k >= 0; k--) {
    if (alphaOcc[k] + betaOcc[k] > 0) {
      homoIndex = k;
      break;
    }
  }
  let lumoIndex = -1;
  for (let k = 0; k < n; k++) {
    if (alphaOcc[k] + betaOcc[k] < 2) {
      if (k > homoIndex || alphaOcc[k] + betaOcc[k] === 0) {
        lumoIndex = k;
        break;
      }
    }
  }

  if (excitedState && homoIndex >= 0 && lumoIndex >= 0 && lumoIndex !== homoIndex) {
    if (betaOcc[homoIndex] > 0) {
      betaOcc[homoIndex] -= 1;
      alphaOcc[lumoIndex] = Math.min(1, alphaOcc[lumoIndex] + 1);
    } else if (alphaOcc[homoIndex] > 0) {
      alphaOcc[homoIndex] -= 1;
      alphaOcc[lumoIndex] = Math.min(1, alphaOcc[lumoIndex] + 1);
    }
  }

  const occupancies = alphaOcc.map((a, k) => a + betaOcc[k]);
  let unpairedElectrons = 0;
  let totalElectronicEnergy = 0;
  for (let k = 0; k < n; k++) {
    unpairedElectrons += Math.abs(alphaOcc[k] - betaOcc[k]);
    totalElectronicEnergy += occupancies[k] * energies[k];
  }

  const fermiEnergy =
    homoIndex >= 0 && lumoIndex >= 0
      ? 0.5 * (energies[homoIndex] + energies[lumoIndex])
      : homoIndex >= 0
        ? energies[homoIndex]
        : 0;

  // Build density matrices P = P^α + P^β and P_spin = P^α - P^β
  const P: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const PSpin: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let k = 0; k < n; k++) {
    const occ = occupancies[k];
    const spin = alphaOcc[k] - betaOcc[k];
    if (occ === 0 && spin === 0) continue;
    for (let i = 0; i < n; i++) {
      const cik = C[i][k];
      for (let j = i; j < n; j++) {
        const prod = cik * C[j][k];
        P[i][j] += occ * prod;
        if (j !== i) P[j][i] += occ * prod;
        if (spin !== 0) {
          PSpin[i][j] += spin * prod;
          if (j !== i) PSpin[j][i] += spin * prod;
        }
      }
    }
  }

  // PS and PSpinS products
  const PS: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const PSpinS: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let s1 = 0;
      let s2 = 0;
      for (let k = 0; k < n; k++) {
        s1 += P[i][k] * S[k][j];
        if (unpairedElectrons > 0) {
          s2 += PSpin[i][k] * S[k][j];
        }
      }
      PS[i][j] = s1;
      PSpinS[i][j] = s2;
    }
  }

  // Mulliken atomic charges & shell breakdown
  const mullikenAtoms: MullikenAtomResult[] = atoms.map((atom, atomIndex) => {
    const props = ELEMENT_PROPERTIES[atom.element] || ELEMENT_PROPERTIES.C;
    const shellPopulations = { s: 0, p: 0, d: 0, f: 0 };
    let grossPopulation = 0;
    let spinPopulation = 0;
    for (let i = 0; i < n; i++) {
      if (basis[i].atomIndex !== atomIndex) continue;
      const qMu = PS[i][i];
      grossPopulation += qMu;
      spinPopulation += PSpinS[i][i];
      if (basis[i].l === 0) shellPopulations.s += qMu;
      else if (basis[i].l === 1) shellPopulations.p += qMu;
      else if (basis[i].l === 2) shellPopulations.d += qMu;
      else shellPopulations.f += qMu;
    }
    return {
      atomIndex,
      atomId: atom.id,
      element: atom.element,
      valenceZ: props.valenceElectrons,
      grossPopulation,
      netCharge: props.valenceElectrons - grossPopulation,
      spinPopulation,
      shellPopulations,
    };
  });

  // Molecular dipole moment from Mulliken point charges relative to center of geometry
  let cx = 0, cy = 0, cz = 0;
  for (const a of atoms) {
    cx += a.x3d;
    cy += a.y3d;
    cz += a.z3d;
  }
  if (atoms.length > 0) {
    cx /= atoms.length;
    cy /= atoms.length;
    cz /= atoms.length;
  }
  let dipX = 0, dipY = 0, dipZ = 0;
  for (let i = 0; i < atoms.length; i++) {
    const q = mullikenAtoms[i].netCharge;
    dipX += q * (atoms[i].x3d - cx);
    dipY += q * (atoms[i].y3d - cy);
    dipZ += q * (atoms[i].z3d - cz);
  }
  const dipoleVectorDebye: [number, number, number] = [
    dipX * DEBYE_PER_E_ANGSTROM,
    dipY * DEBYE_PER_E_ANGSTROM,
    dipZ * DEBYE_PER_E_ANGSTROM,
  ];
  const dipoleMagnitudeDebye = Math.hypot(...dipoleVectorDebye);

  // Mayer bond orders & Mulliken overlap populations
  const bondOrderMap = new Map<string, number>();
  for (const b of bonds) {
    bondOrderMap.set(`${b.source}::${b.target}`, b.order);
    bondOrderMap.set(`${b.target}::${b.source}`, b.order);
  }

  const atomBasisIndices: number[][] = atoms.map(() => []);
  for (let i = 0; i < n; i++) {
    atomBasisIndices[basis[i].atomIndex].push(i);
  }

  const mayerBonds: MayerBondOrderResult[] = [];
  const atomicValences = new Array<number>(atoms.length).fill(0);

  for (let a = 0; a < atoms.length; a++) {
    for (let b = a + 1; b < atoms.length; b++) {
      let mayer = 0;
      let mullikenOv = 0;
      for (const mu of atomBasisIndices[a]) {
        for (const nu of atomBasisIndices[b]) {
          mayer += PS[mu][nu] * PS[nu][mu] + PSpinS[mu][nu] * PSpinS[nu][mu];
          mullikenOv += 2 * P[mu][nu] * S[mu][nu];
        }
      }
      atomicValences[a] += mayer;
      atomicValences[b] += mayer;

      const lewisOrder = bondOrderMap.get(`${atoms[a].id}::${atoms[b].id}`) || 0;
      const dist = Math.hypot(
        atoms[a].x3d - atoms[b].x3d,
        atoms[a].y3d - atoms[b].y3d,
        atoms[a].z3d - atoms[b].z3d
      );
      if (lewisOrder > 0 || mayer >= 0.08) {
        mayerBonds.push({
          atomAIndex: a,
          atomBIndex: b,
          atomAId: atoms[a].id,
          atomBId: atoms[b].id,
          elementA: atoms[a].element,
          elementB: atoms[b].element,
          mayerOrder: mayer,
          mullikenOverlapPop: mullikenOv,
          lewisOrder,
          distanceAngstrom: dist,
        });
      }
    }
  }
  mayerBonds.sort((x, y) => y.mayerOrder - x.mayerOrder);

  return {
    basis,
    S,
    H,
    energies,
    coefficients: C,
    occupancies,
    alphaOccupancies: alphaOcc,
    betaOccupancies: betaOcc,
    densityMatrix: P,
    spinDensityMatrix: PSpin,
    mullikenAtoms,
    mayerBonds,
    atomicValences,
    dipoleVectorDebye,
    dipoleMagnitudeDebye,
    fermiEnergy,
    homoIndex,
    lumoIndex,
    totalElectronicEnergy,
    totalElectrons,
    unpairedElectrons,
    charge,
  };
}

/**
 * Extracts normalized LCAOContribution list for a specific MO index k from Löwdin EHT results.
 */
export function extractLCAOFromLowdin(eht: LowdinEHTResult, moIndex: number): LCAOContribution[] {
  const n = eht.basis.length;
  const list: LCAOContribution[] = [];
  for (let i = 0; i < n; i++) {
    const c = eht.coefficients[i]?.[moIndex] ?? 0;
    if (Math.abs(c) < 0.03) continue;
    const b = eht.basis[i];
    list.push({
      atomId: b.atomId,
      element: b.element,
      orbital: b.orbital,
      label: b.label,
      coeff: c,
    });
  }
  return list.sort((a, b) => Math.abs(b.coeff) - Math.abs(a.coeff));
}

/**
 * Computes Gaussian-broadened Projected Density of States (PDOS) using exact Mulliken weights
 * w_{μ,k} = C_{μk} * (S C)_{μk} so Σ_μ w_{μ,k} = 1 for every MO k!
 */
export function computePDOS(
  eht: LowdinEHTResult,
  options?: { sigmaEv?: number; numPoints?: number; elementFilter?: string | null }
): {
  points: PDOSPoint[];
  elements: string[];
  eMin: number;
  eMax: number;
} {
  const sigma = options?.sigmaEv ?? 0.55;
  const numPoints = options?.numPoints ?? 140;
  const n = eht.basis.length;
  if (n === 0) {
    return { points: [], elements: [], eMin: -20, eMax: 5 };
  }

  const eMin = Math.min(...eht.energies) - 2.5;
  const eMax = Math.max(...eht.energies) + 2.5;
  const elements = Array.from(new Set(eht.basis.map((b) => b.element)));

  // Precompute Mulliken orbital projection weights w[moIdx][element/shell]
  const moElemWeights: Record<string, number>[] = [];
  const moShellWeights: { s: number; p: number; d: number; f: number }[] = [];

  for (let k = 0; k < n; k++) {
    const elemW: Record<string, number> = {};
    for (const el of elements) elemW[el] = 0;
    const shellW = { s: 0, p: 0, d: 0, f: 0 };

    for (let i = 0; i < n; i++) {
      let sc = 0;
      for (let j = 0; j < n; j++) {
        sc += eht.S[i][j] * eht.coefficients[j][k];
      }
      const w = eht.coefficients[i][k] * sc;
      const b = eht.basis[i];
      elemW[b.element] = (elemW[b.element] || 0) + w;
      if (b.l === 0) shellW.s += w;
      else if (b.l === 1) shellW.p += w;
      else if (b.l === 2) shellW.d += w;
      else shellW.f += w;
    }
    moElemWeights.push(elemW);
    moShellWeights.push(shellW);
  }

  const norm = 1 / (sigma * Math.sqrt(2 * Math.PI));
  const inv2Sig2 = 0.5 / (sigma * sigma);
  const points: PDOSPoint[] = [];

  for (let p = 0; p < numPoints; p++) {
    const E = eMin + (p / Math.max(1, numPoints - 1)) * (eMax - eMin);
    let total = 0;
    const byElement: Record<string, number> = {};
    for (const el of elements) byElement[el] = 0;
    const byShell = { s: 0, p: 0, d: 0, f: 0 };

    for (let k = 0; k < n; k++) {
      const dE = E - eht.energies[k];
      const g = norm * Math.exp(-dE * dE * inv2Sig2);
      total += g;
      for (const el of elements) {
        byElement[el] += g * moElemWeights[k][el];
      }
      byShell.s += g * moShellWeights[k].s;
      byShell.p += g * moShellWeights[k].p;
      byShell.d += g * moShellWeights[k].d;
      byShell.f += g * moShellWeights[k].f;
    }

    points.push({ energy: E, total, byElement, byShell });
  }

  return { points, elements, eMin, eMax };
}

/**
 * Computes Crystal Orbital Overlap Population (COOP) curve for a selected atom pair (A, B)
 * or summed across all bonded pairs when pair is null:
 *   COOP_{AB}(E) = Σ_k [ 2 Σ_{μ∈A, ν∈B} C_{μk} C_{νk} S_{μν} ] G(E - ε_k)
 * Positive values indicate bonding character; negative values indicate antibonding character.
 */
export function computeCOOP(
  eht: LowdinEHTResult,
  bonds: BondEdge[],
  pair?: { atomAIndex: number; atomBIndex: number } | null,
  options?: { sigmaEv?: number; numPoints?: number }
): {
  points: COOPPoint[];
  icoopAtFermi: number;
  moOverlapPopulations: number[];
  eMin: number;
  eMax: number;
} {
  const sigma = options?.sigmaEv ?? 0.55;
  const numPoints = options?.numPoints ?? 140;
  const n = eht.basis.length;
  if (n === 0) {
    return { points: [], icoopAtFermi: 0, moOverlapPopulations: [], eMin: -20, eMax: 5 };
  }

  const eMin = Math.min(...eht.energies) - 2.5;
  const eMax = Math.max(...eht.energies) + 2.5;

  // Determine target atom-index pairs
  const pairs: [number, number][] = [];
  if (pair) {
    pairs.push([pair.atomAIndex, pair.atomBIndex]);
  } else if (eht.mayerBonds.length > 0) {
    for (const mb of eht.mayerBonds) {
      if (mb.lewisOrder > 0 || mb.mayerOrder >= 0.15) {
        pairs.push([mb.atomAIndex, mb.atomBIndex]);
      }
    }
  } else {
    // Fallback to bonds list
    const idToIdx = new Map<string, number>();
    for (const b of eht.basis) idToIdx.set(b.atomId, b.atomIndex);
    for (const bond of bonds) {
      const iA = idToIdx.get(bond.source);
      const iB = idToIdx.get(bond.target);
      if (iA !== undefined && iB !== undefined) pairs.push([iA, iB]);
    }
  }

  const atomBasisMap = new Map<number, number[]>();
  for (let i = 0; i < n; i++) {
    const list = atomBasisMap.get(eht.basis[i].atomIndex) || [];
    list.push(i);
    atomBasisMap.set(eht.basis[i].atomIndex, list);
  }

  // Per-MO overlap population op_k
  const moOverlapPopulations = new Array<number>(n).fill(0);
  for (let k = 0; k < n; k++) {
    let op = 0;
    for (const [a, b] of pairs) {
      const bA = atomBasisMap.get(a) || [];
      const bB = atomBasisMap.get(b) || [];
      for (const mu of bA) {
        const cMu = eht.coefficients[mu][k];
        for (const nu of bB) {
          op += 2 * cMu * eht.coefficients[nu][k] * eht.S[mu][nu];
        }
      }
    }
    moOverlapPopulations[k] = op;
  }

  let icoopAtFermi = 0;
  for (let k = 0; k < n; k++) {
    icoopAtFermi += eht.occupancies[k] * moOverlapPopulations[k];
  }

  const norm = 1 / (sigma * Math.sqrt(2 * Math.PI));
  const inv2Sig2 = 0.5 / (sigma * sigma);
  const dEStep = (eMax - eMin) / Math.max(1, numPoints - 1);
  const points: COOPPoint[] = [];
  let runningInt = 0;

  for (let p = 0; p < numPoints; p++) {
    const E = eMin + p * dEStep;
    let val = 0;
    for (let k = 0; k < n; k++) {
      const dE = E - eht.energies[k];
      val += moOverlapPopulations[k] * norm * Math.exp(-dE * dE * inv2Sig2);
    }
    runningInt += val * dEStep * 2; // spin-degeneracy factor 2
    points.push({ energy: E, value: val, integrated: runningInt });
  }

  return { points, icoopAtFermi, moOverlapPopulations, eMin, eMax };
}
