import { AtomNode, LCAOContribution, MolecularOrbital, OrbitalBasisType } from './moleculesData';
import {
  ANGSTROM_TO_BOHR,
  LowdinEHTResult,
  getSlaterExponentForBasis,
  symmetricJacobiEigen,
} from './quantumEngineV3';
import { GaussianCubeData } from './fileIO';
import { IsosurfaceMeshData } from './moIsosurface';

export type VolumetricFieldMode = 'mo' | 'density' | 'mep' | 'spin' | 'cube';
export type SliceOrientation = 'xy' | 'xz' | 'yz';

export interface ColoredIsosurfaceMeshData extends IsosurfaceMeshData {
  colors?: Float32Array; // Optional per-vertex RGB [0..1] (used for MEP-on-density)
}

export interface VolumetricSurfaceResult {
  positive: ColoredIsosurfaceMeshData;
  negative: ColoredIsosurfaceMeshData;
  mepRange?: { min: number; max: number };
  hasData: boolean;
}

export interface Slice2DResult {
  width: number;
  height: number;
  values: Float32Array;
  rgba: Uint8ClampedArray;
  minVal: number;
  maxVal: number;
  absMax: number;
  spanU: number;
  spanV: number;
  center: [number, number, number];
  orientation: SliceOrientation;
  offset: number;
}

// ============================================================================
// Real Normalized Slater-Type Orbital (STO) Spatial Evaluation
// ============================================================================

const FACT_SMALL = [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600];

/**
 * Evaluates a real normalized Slater-type orbital φ_{n,l,m}(r) in atomic units (a0^{-3/2})
 * at Cartesian displacement (dx, dy, dz) in Angstroms from the atomic center.
 */
export function evalRealSTOBasis(
  element: AtomNode['element'],
  orbital: OrbitalBasisType,
  dxAng: number,
  dyAng: number,
  dzAng: number
): number {
  const { n, zeta } = getSlaterExponentForBasis(element, orbital);
  const dx = dxAng * ANGSTROM_TO_BOHR;
  const dy = dyAng * ANGSTROM_TO_BOHR;
  const dz = dzAng * ANGSTROM_TO_BOHR;
  const r2 = dx * dx + dy * dy + dz * dz;
  const r = Math.sqrt(r2);
  if (r > 14.0) return 0;

  const twoZ = 2 * zeta;
  const radialNorm = Math.pow(twoZ, n + 0.5) / Math.sqrt(FACT_SMALL[2 * n] || 479001600);
  const expPart = Math.exp(-zeta * r);

  // s-orbital (l = 0): Y_00 = 1 / sqrt(4π) = 0.28209479
  if (orbital === 's') {
    const rPow = n === 1 ? 1 : Math.pow(r, n - 1);
    return radialNorm * rPow * expPart * 0.28209479177387814;
  }

  // p-orbitals (l = 1): Y_1m = sqrt(3/(4π)) * (x|y|z)/r = 0.48860251 * coord/r
  if (orbital === 'px' || orbital === 'py' || orbital === 'pz') {
    const rPow = n <= 2 ? 1 : Math.pow(r, n - 2);
    const coord = orbital === 'px' ? dx : orbital === 'py' ? dy : dz;
    return radialNorm * rPow * expPart * 0.4886025119029199 * coord;
  }

  // d-orbitals (l = 2)
  if (orbital.startsWith('d')) {
    const rPow = n <= 3 ? 1 : Math.pow(r, n - 3);
    let ang = 0;
    if (orbital === 'dz2') ang = 0.31539156525252005 * (2 * dz * dz - dx * dx - dy * dy);
    else if (orbital === 'dx2y2') ang = 0.5462742152960396 * (dx * dx - dy * dy);
    else if (orbital === 'dxy') ang = 1.0925484305920792 * dx * dy;
    else if (orbital === 'dxz') ang = 1.0925484305920792 * dx * dz;
    else if (orbital === 'dyz') ang = 1.0925484305920792 * dy * dz;
    return radialNorm * rPow * expPart * ang;
  }

  // f-orbitals (l = 3)
  const rPow = n <= 4 ? 1 : Math.pow(r, n - 4);
  let ang = 0;
  if (orbital === 'fz3') ang = 0.3731763325901154 * dz * (2 * dz * dz - 3 * dx * dx - 3 * dy * dy);
  else if (orbital === 'fxz2') ang = 0.4570457994644658 * dx * (4 * dz * dz - dx * dx - dy * dy);
  else if (orbital === 'fyz2') ang = 0.4570457994644658 * dy * (4 * dz * dz - dx * dx - dy * dy);
  else if (orbital === 'fz_x2y2') ang = 1.445305721320277 * dz * (dx * dx - dy * dy);
  else if (orbital === 'fxyz') ang = 2.890611442640554 * dx * dy * dz;
  else if (orbital === 'fx_x2_3y2') ang = 0.5900435899266435 * dx * (dx * dx - 3 * dy * dy);
  else if (orbital === 'fy_3x2_y2') ang = 0.5900435899266435 * dy * (3 * dx * dx - dy * dy);
  return radialNorm * rPow * expPart * ang;
}

/**
 * Evaluates an MO wavefunction ψ(r) using real shell-specific Slater-type orbitals.
 */
export function evalRealSTOWavefunction(
  lcao: LCAOContribution[],
  atoms: AtomNode[],
  x: number,
  y: number,
  z: number
): number {
  const atomById = new Map(atoms.map((a) => [a.id, a]));
  let psi = 0;
  for (const term of lcao) {
    if (Math.abs(term.coeff) < 0.02) continue;
    const atom = atomById.get(term.atomId);
    if (!atom) continue;
    psi += term.coeff * evalRealSTOBasis(atom.element, term.orbital, x - atom.x3d, y - atom.y3d, z - atom.z3d);
  }
  return psi;
}

/**
 * Evaluates total valence electron density ρ(r) = Σ_k n_k |ψ_k(r)|^2 and spin density ρ_α - ρ_β
 * using the Löwdin generalized EHT density matrix & real STO basis.
 */
export function evalDensityAndSpinAtPoint(
  eht: LowdinEHTResult,
  atoms: AtomNode[],
  x: number,
  y: number,
  z: number
): { density: number; spinDensity: number } {
  const n = eht.basis.length;
  const phi = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const b = eht.basis[i];
    const a = atoms[b.atomIndex];
    const dx = x - a.x3d;
    const dy = y - a.y3d;
    const dz = z - a.z3d;
    if (dx * dx + dy * dy + dz * dz > 28) continue;
    phi[i] = evalRealSTOBasis(b.element, b.orbital, dx, dy, dz);
  }

  let density = 0;
  let spinDensity = 0;
  const hasSpin = eht.unpairedElectrons > 0;
  for (let i = 0; i < n; i++) {
    const pi = phi[i];
    if (pi === 0) continue;
    // Diagonal
    density += eht.densityMatrix[i][i] * pi * pi;
    if (hasSpin) spinDensity += eht.spinDensityMatrix[i][i] * pi * pi;
    // Off-diagonal
    for (let j = i + 1; j < n; j++) {
      const pj = phi[j];
      if (pj === 0) continue;
      const prod = 2 * pi * pj;
      density += eht.densityMatrix[i][j] * prod;
      if (hasSpin) spinDensity += eht.spinDensityMatrix[i][j] * prod;
    }
  }
  return { density: Math.max(0, density), spinDensity };
}

// ============================================================================
// Quantum Distributed Multipole Analysis (SCC-EHT + CM5 + Analytical STO DMA)
// ============================================================================

const HUBBARD_U_EV: Partial<Record<AtomNode['element'], number>> = {
  H: 12.85, He: 24.0,
  Li: 4.8, Be: 8.8, B: 8.0, C: 10.0, N: 11.7, O: 13.4, F: 14.9, Ne: 21.0,
  Na: 4.6, Mg: 7.6, Al: 5.5, Si: 7.0, P: 8.2, S: 9.0, Cl: 9.9, Ar: 15.0,
  K: 4.3, Ca: 6.1, Ti: 6.8, V: 6.7, Cr: 6.7, Mn: 7.4, Fe: 7.9, Co: 7.9,
  Ni: 7.6, Cu: 7.7, Zn: 9.4, Ga: 5.8, Ge: 6.8, As: 7.9, Se: 8.2, Br: 8.8,
  Ru: 7.3, Rh: 7.5, Pd: 8.3, Ag: 7.6, Sn: 6.3, Sb: 7.3, Te: 7.5, I: 7.6,
  Xe: 12.1, Ir: 8.0, Pt: 8.2, Au: 8.5,
};

const PAULING_EN_TABLE: Partial<Record<AtomNode['element'], number>> = {
  H: 2.20, He: 2.50,
  Li: 0.98, Be: 1.57, B: 2.04, C: 2.55, N: 3.04, O: 3.44, F: 3.98, Ne: 3.00,
  Na: 0.93, Mg: 1.31, Al: 1.61, Si: 1.90, P: 2.19, S: 2.58, Cl: 3.16, Ar: 2.80,
  K: 0.82, Ca: 1.00, Ti: 1.54, V: 1.63, Cr: 1.66, Mn: 1.55, Fe: 1.83, Co: 1.88,
  Ni: 1.91, Cu: 1.90, Zn: 1.65, Ga: 1.81, Ge: 2.01, As: 2.18, Se: 2.55, Br: 2.96,
  Ru: 2.20, Rh: 2.28, Pd: 2.20, Ag: 1.93, Sn: 1.96, Sb: 2.05, Te: 2.10, I: 2.66,
  Xe: 2.60, Ce: 1.12, Ir: 2.20, Pt: 2.28, Au: 2.54,
};

const FACT_DMA = [
  1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
  6227020800, 87178291200, 1307674368000,
];

function factDMA(n: number): number {
  if (n <= 1) return 1;
  if (n < FACT_DMA.length) return FACT_DMA[n];
  let r = FACT_DMA[FACT_DMA.length - 1];
  for (let i = FACT_DMA.length; i <= n; i++) r *= i;
  return r;
}

function expPolySum(m: number, t: number): number {
  let term = 1;
  let sum = 1;
  for (let j = 1; j <= m; j++) {
    term *= t / j;
    sum += term;
  }
  return sum;
}

/**
 * Exact closed-form radial Slater-type orbital multipole potential integral:
 *   I_L(r) = N1 N2 [ r^{-(L+1)} ∫_0^r r'^{k+L} e^{-αr'} dr' + r^L ∫_r^∞ r'^{k-L-1} e^{-αr'} dr' ]
 * where k = n1 + n2, α = ζ1 + ζ2.
 */
function evalSTORadialMultipolePrecomp(
  r: number,
  L: 1 | 2,
  n1n2Norm: number,
  alpha: number,
  k: number
): number {
  if (r < 1e-5) return 0;
  const mLow = k + L;
  const lowCoeff = factDMA(mLow) / Math.pow(alpha, mLow + 1);
  const t = alpha * r;
  if (t > 36) {
    return n1n2Norm * (lowCoeff / (L === 1 ? r * r : r * r * r));
  }
  const emt = Math.exp(-t);
  const lowInt = lowCoeff * (1 - emt * expPolySum(mLow, t));
  const mHigh = k - L - 1;
  const highInt =
    mHigh >= 0
      ? (factDMA(mHigh) / Math.pow(alpha, mHigh + 1)) * emt * expPolySum(mHigh, t)
      : 0;
  const rPowLow = L === 1 ? r * r : r * r * r;
  const rPowHigh = L === 1 ? r : r * r;
  return n1n2Norm * (lowInt / rPowLow + rPowHigh * highInt);
}

interface MolecularDMAAtom {
  x3d: number;
  y3d: number;
  z3d: number;
  q: number; // Self-consistent CM5 monopole charge (e)
  // L=1 s-p lone-pair hybridization dipole
  hasSP: boolean;
  spNorm: number;
  spAlpha: number;
  spK: number;
  spCoeffX: number; // pre-multiplied by nScale * (2 / sqrt(3)) * P_{s, px}
  spCoeffY: number;
  spCoeffZ: number;
  // L=2 p-p valence quadrupole (captures π-clouds and halogen/chalcogen σ-holes)
  hasP: boolean;
  pNorm: number;
  pAlpha: number;
  pK: number;
  pXX: number;
  pYY: number;
  pZZ: number;
  pXY: number;
  pXZ: number;
  pYZ: number;
  pTrace: number;
}

interface MolecularDMA {
  atoms: MolecularDMAAtom[];
}

const molecularDMACache = new WeakMap<LowdinEHTResult, MolecularDMA>();

/**
 * Builds a cached Distributed Multipole Analysis (DMA) model for a molecule from its
 * Löwdin generalized EHT wavefunction:
 *   1. Self-Consistent Charge EHT (SCC-EHT) with Hubbard-U Ohno-Klopman Coulomb feedback
 *      removes raw non-iterative EHT charge exaggeration (e.g., over-polarized O/N/F and
 *      artificially positive Period-3 S in thiophene).
 *   2. Universal charge-conserving CM5 bond-dipole correction fixes minimal-basis X-H
 *      electronegativity bias while strictly preserving Σ_A q_A = Q_mol.
 *   3. Exact analytical 1-center STO L=1 (s-p lone-pair dipole) and L=2 (p-p π-cloud &
 *      σ-hole quadrupole) tensors capture anisotropic electrostatic topology missed by
 *      spherical point charges (aromatic π-rings, N2 quadrupole, Cl σ-hole, lone pairs).
 */
function getOrBuildMolecularDMA(eht: LowdinEHTResult, atoms: AtomNode[]): MolecularDMA {
  const cached = molecularDMACache.get(eht);
  if (cached && cached.atoms.length === atoms.length) {
    return cached;
  }

  const n = eht.basis.length;
  const N = atoms.length;
  if (n === 0 || N === 0) {
    const empty: MolecularDMA = { atoms: [] };
    molecularDMACache.set(eht, empty);
    return empty;
  }

  const S = eht.S;
  const basis = eht.basis;
  const occ = eht.occupancies;

  // 1. Löwdin orthogonalizer X = S^{-1/2}
  const sEig = symmetricJacobiEigen(S);
  const invSqrt = sEig.eigenvalues.map((v) => 1.0 / Math.sqrt(Math.max(1e-5, v)));
  const U = sEig.eigenvectors;
  const X: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      let sx = 0;
      for (let k = 0; k < n; k++) sx += U[i][k] * invSqrt[k] * U[j][k];
      X[i][j] = sx;
      X[j][i] = sx;
    }
  }

  // 2. Ohno-Klopman Coulomb matrix γ_AB (in eV) for SCC-EHT charge equilibration
  const KE = 14.3996;
  const gamma: number[][] = Array.from({ length: N }, () => new Array<number>(N).fill(0));
  for (let a = 0; a < N; a++) {
    const ua = HUBBARD_U_EV[atoms[a].element] ?? 9.2;
    gamma[a][a] = ua;
    for (let b = a + 1; b < N; b++) {
      const ub = HUBBARD_U_EV[atoms[b].element] ?? 9.2;
      const r = Math.hypot(
        atoms[a].x3d - atoms[b].x3d,
        atoms[a].y3d - atoms[b].y3d,
        atoms[a].z3d - atoms[b].z3d
      );
      const a0 = KE / (0.5 * (ua + ub));
      const g = KE / Math.sqrt(r * r + a0 * a0);
      gamma[a][b] = g;
      gamma[b][a] = g;
    }
  }

  // Period-3+ chalcogen/pnictogen valence s-orbital VSIE calibration prevents artificial
  // 3s-2s electron loss to carbon while leaving H 1s and Period-2 EHT parameters untouched.
  const hii0 = basis.map((b) => {
    if (b.element === 'S' && b.l === 0) return -22.2;
    if (b.element === 'S' && b.l === 1) return -13.7;
    if (b.element === 'P' && b.l === 0) return -19.8;
    if (b.element === 'Se' && b.l === 0) return -21.2;
    return b.hii;
  });

  const H0: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    H0[i][i] = hii0[i];
    for (let j = i + 1; j < n; j++) {
      const hij = 0.5 * 1.75 * (hii0[i] + hii0[j]) * S[i][j];
      H0[i][j] = hij;
      H0[j][i] = hij;
    }
  }

  // Pre-group basis indices per atom
  const atomBasisIndices: number[][] = Array.from({ length: N }, () => []);
  for (let i = 0; i < n; i++) {
    atomBasisIndices[basis[i].atomIndex].push(i);
  }

  // 3. Damped Self-Consistent Charge (SCC) loop
  const qAtom = new Array<number>(N).fill(0);
  const vShift = new Array<number>(N).fill(0);
  const H: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const HX: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  const HPrime: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  let lastCPrime: number[][] = [];
  const maxIter = n > 68 ? 5 : 12;
  const mix = n > 68 ? 0.48 : 0.38;

  for (let iter = 0; iter < maxIter; iter++) {
    for (let a = 0; a < N; a++) {
      let s = 0;
      for (let b = 0; b < N; b++) s += gamma[a][b] * qAtom[b];
      vShift[a] = -s;
    }
    for (let i = 0; i < n; i++) {
      const vi = vShift[basis[i].atomIndex];
      H[i][i] = H0[i][i] + vi;
      for (let j = i + 1; j < n; j++) {
        const vj = vShift[basis[j].atomIndex];
        const hij = H0[i][j] + 0.5 * S[i][j] * (vi + vj);
        H[i][j] = hij;
        H[j][i] = hij;
      }
    }
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let k = 0; k < n; k++) s += H[i][k] * X[k][j];
        HX[i][j] = s;
      }
    }
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        let s = 0;
        for (let k = 0; k < n; k++) s += X[i][k] * HX[k][j];
        HPrime[i][j] = s;
        HPrime[j][i] = s;
      }
    }

    const eig = symmetricJacobiEigen(HPrime);
    lastCPrime = eig.eigenvectors;

    let maxDiff = 0;
    for (let a = 0; a < N; a++) {
      let pop = 0;
      for (const i of atomBasisIndices[a]) {
        for (let k = 0; k < n; k++) {
          if (occ[k] > 0) pop += occ[k] * lastCPrime[i][k] * lastCPrime[i][k];
        }
      }
      const zVal = eht.mullikenAtoms[a]?.valenceZ ?? 0;
      const qNew = zVal - pop;
      const diff = Math.abs(qNew - qAtom[a]);
      if (diff > maxDiff) maxDiff = diff;
      qAtom[a] = (1 - mix) * qAtom[a] + mix * qNew;
    }
    if (maxDiff < 5e-4) break;
  }

  // Back-transform C = X C' once to obtain true AO density matrix P for 1-center multipoles
  const C: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < n; k++) {
      let s = 0;
      for (let m = 0; m < n; m++) s += X[i][m] * lastCPrime[m][k];
      C[i][k] = s;
    }
  }

  // 4. Universal charge-conserving CM5 bond-dipole correction (exact Σ q_A = Q_mol)
  const qCM5 = qAtom.slice();
  for (const b of eht.mayerBonds) {
    const enA = PAULING_EN_TABLE[b.elementA] ?? 2.2;
    const enB = PAULING_EN_TABLE[b.elementB] ?? 2.2;
    const bo = Math.min(2.5, Math.max(0, b.mayerOrder));
    if (bo < 1e-4) continue;
    let dq = 0;
    if (b.elementA === 'H' || b.elementB === 'H') {
      const hIsA = b.elementA === 'H';
      const otherEN = hIsA ? enB : enA;
      const shiftH =
        (otherEN >= 2.2
          ? 0.125 * Math.tanh((otherEN - 2.2) / 0.25)
          : -0.1 * Math.tanh((2.2 - otherEN) / 0.25)) * bo;
      dq = hIsA ? shiftH : -shiftH;
    } else {
      dq = 0.025 * (enB - enA) * bo;
    }
    qCM5[b.atomAIndex] += dq;
    qCM5[b.atomBIndex] -= dq;
  }

  // Helper to compute single P[i][j] element on demand (only needed for 1-center s/p blocks!)
  const evalPij = (i: number, j: number): number => {
    let s = 0;
    for (let k = 0; k < n; k++) {
      if (occ[k] > 0) s += occ[k] * C[i][k] * C[j][k];
    }
    return s;
  };

  // 5. Assemble precomputed per-atom DMA records
  const TWO_OVER_SQRT3 = 2 / Math.sqrt(3);
  const dmaAtoms: MolecularDMAAtom[] = atoms.map((atom, a) => {
    const idxs = atomBasisIndices[a];
    let sIdx = -1, pxIdx = -1, pyIdx = -1, pzIdx = -1;
    for (const idx of idxs) {
      const orb = basis[idx].orbital;
      if (orb === 's') sIdx = idx;
      else if (orb === 'px') pxIdx = idx;
      else if (orb === 'py') pyIdx = idx;
      else if (orb === 'pz') pzIdx = idx;
    }

    const hasP = pxIdx >= 0 && pyIdx >= 0 && pzIdx >= 0;
    const hasSP = sIdx >= 0 && hasP;

    let spNorm = 0, spAlpha = 1, spK = 4, spCoeffX = 0, spCoeffY = 0, spCoeffZ = 0;
    if (hasSP) {
      const bs = basis[sIdx];
      const bp = basis[pxIdx];
      const n1 = Math.pow(2 * bs.zeta, bs.n + 0.5) / Math.sqrt(factDMA(2 * bs.n));
      const n2 = Math.pow(2 * bp.zeta, bp.n + 0.5) / Math.sqrt(factDMA(2 * bp.n));
      spNorm = n1 * n2;
      spAlpha = bs.zeta + bp.zeta;
      spK = bs.n + bp.n;
      // Radial-node damping for nodeless STOs: real 2s has 1 node (~0.55x), 3s+ has 2+ nodes (~0.22x)
      const nScale = bs.n >= 3 ? 0.22 : 0.55;
      const factor = nScale * TWO_OVER_SQRT3;
      spCoeffX = factor * evalPij(sIdx, pxIdx);
      spCoeffY = factor * evalPij(sIdx, pyIdx);
      spCoeffZ = factor * evalPij(sIdx, pzIdx);
    }

    let pNorm = 0, pAlpha = 1, pK = 4;
    let pXX = 0, pYY = 0, pZZ = 0, pXY = 0, pXZ = 0, pYZ = 0, pTrace = 0;
    if (hasP) {
      const bp = basis[pxIdx];
      const np = Math.pow(2 * bp.zeta, bp.n + 0.5) / Math.sqrt(factDMA(2 * bp.n));
      pNorm = np * np;
      pAlpha = 2 * bp.zeta;
      pK = 2 * bp.n;
      pXX = evalPij(pxIdx, pxIdx);
      pYY = evalPij(pyIdx, pyIdx);
      pZZ = evalPij(pzIdx, pzIdx);
      pXY = evalPij(pxIdx, pyIdx);
      pXZ = evalPij(pxIdx, pzIdx);
      pYZ = evalPij(pyIdx, pzIdx);
      pTrace = pXX + pYY + pZZ;
    }

    return {
      x3d: atom.x3d,
      y3d: atom.y3d,
      z3d: atom.z3d,
      q: qCM5[a],
      hasSP,
      spNorm,
      spAlpha,
      spK,
      spCoeffX,
      spCoeffY,
      spCoeffZ,
      hasP,
      pNorm,
      pAlpha,
      pK,
      pXX,
      pYY,
      pZZ,
      pXY,
      pXZ,
      pYZ,
      pTrace,
    };
  });

  const dma: MolecularDMA = { atoms: dmaAtoms };
  molecularDMACache.set(eht, dma);
  return dma;
}

/**
 * Evaluates the Molecular Electrostatic Potential V_ESP(r) (in atomic units Hartree/e)
 * via Quantum Distributed Multipole Analysis (SCC-EHT + CM5 monopoles + analytical STO
 * L=1 s-p lone-pair dipoles + L=2 p-p π-cloud / σ-hole quadrupoles).
 * Negative V(r) = electron-rich nucleophilic region (red); Positive V(r) = electrophilic (blue).
 */
export function evalElectrostaticPotentialAtPoint(
  eht: LowdinEHTResult,
  atoms: AtomNode[],
  x: number,
  y: number,
  z: number
): number {
  const dma = getOrBuildMolecularDMA(eht, atoms);
  const dmaAtoms = dma.atoms;
  let v = 0;

  for (let i = 0; i < dmaAtoms.length; i++) {
    const a = dmaAtoms[i];
    const dx = (x - a.x3d) * ANGSTROM_TO_BOHR;
    const dy = (y - a.y3d) * ANGSTROM_TO_BOHR;
    const dz = (z - a.z3d) * ANGSTROM_TO_BOHR;
    const r2 = dx * dx + dy * dy + dz * dz;
    const r = Math.sqrt(r2);

    // L = 0: Self-consistent CM5 monopole with gentle nuclear softening
    v += a.q / Math.sqrt(r2 + 0.2);

    if (r < 1e-5) continue;
    const invR = 1 / r;
    const ux = dx * invR;
    const uy = dy * invR;
    const uz = dz * invR;

    // L = 1: Analytical 1-center s-p hybridization lone-pair dipole
    if (a.hasSP) {
      const pDotU = a.spCoeffX * ux + a.spCoeffY * uy + a.spCoeffZ * uz;
      if (Math.abs(pDotU) > 1e-7) {
        const I1 = evalSTORadialMultipolePrecomp(r, 1, a.spNorm, a.spAlpha, a.spK);
        v -= pDotU * I1;
      }
    }

    // L = 2: Analytical 1-center p-p atomic quadrupole (aromatic π-clouds, σ-holes, N2)
    if (a.hasP) {
      const uPu =
        a.pXX * ux * ux +
        a.pYY * uy * uy +
        a.pZZ * uz * uz +
        2 * (a.pXY * ux * uy + a.pXZ * ux * uz + a.pYZ * uy * uz);
      const quad = 3 * uPu - a.pTrace;
      if (Math.abs(quad) > 1e-7) {
        const I2 = evalSTORadialMultipolePrecomp(r, 2, a.pNorm, a.pAlpha, a.pK);
        v -= 0.2 * quad * I2;
      }
    }
  }

  return v;
}

/**
 * Evaluates an imported Gaussian .cube grid via trilinear interpolation at world point (x, y, z) in Angstroms.
 */
export function sampleGaussianCubeAtPoint(
  cube: GaussianCubeData,
  x: number,
  y: number,
  z: number
): number {
  const dx = cube.vxAngstrom[0] || 0.2;
  const dy = cube.vyAngstrom[1] || 0.2;
  const dz = cube.vzAngstrom[2] || 0.2;
  const fx = (x - cube.originAngstrom[0]) / dx;
  const fy = (y - cube.originAngstrom[1]) / dy;
  const fz = (z - cube.originAngstrom[2]) / dz;
  if (fx < 0 || fy < 0 || fz < 0 || fx >= cube.nx - 1 || fy >= cube.ny - 1 || fz >= cube.nz - 1) {
    return 0;
  }
  const ix = Math.floor(fx);
  const iy = Math.floor(fy);
  const iz = Math.floor(fz);
  const tx = fx - ix;
  const ty = fy - iy;
  const tz = fz - iz;
  const idx = (a: number, b: number, c: number) => (a * cube.ny + b) * cube.nz + c;

  const c000 = cube.values[idx(ix, iy, iz)];
  const c100 = cube.values[idx(ix + 1, iy, iz)];
  const c010 = cube.values[idx(ix, iy + 1, iz)];
  const c110 = cube.values[idx(ix + 1, iy + 1, iz)];
  const c001 = cube.values[idx(ix, iy, iz + 1)];
  const c101 = cube.values[idx(ix + 1, iy, iz + 1)];
  const c011 = cube.values[idx(ix, iy + 1, iz + 1)];
  const c111 = cube.values[idx(ix + 1, iy + 1, iz + 1)];

  const c00 = c000 * (1 - tx) + c100 * tx;
  const c10 = c010 * (1 - tx) + c110 * tx;
  const c01 = c001 * (1 - tx) + c101 * tx;
  const c11 = c011 * (1 - tx) + c111 * tx;
  const c0 = c00 * (1 - ty) + c10 * ty;
  const c1 = c01 * (1 - ty) + c11 * ty;
  return c0 * (1 - tz) + c1 * tz;
}

// ============================================================================
// Marching Tetrahedra Isosurface Extraction for Generic Volumetric Fields
// ============================================================================

const KUHN_TETS: readonly [number, number, number, number][] = [
  [0, 5, 1, 6],
  [0, 1, 2, 6],
  [0, 2, 3, 6],
  [0, 3, 7, 6],
  [0, 7, 4, 6],
  [0, 4, 5, 6],
];

const TET_EDGES: readonly (readonly number[])[] = [
  [],
  [0, 1, 0, 2, 0, 3],
  [1, 0, 1, 3, 1, 2],
  [0, 2, 0, 3, 1, 3, 1, 2],
  [2, 0, 2, 1, 2, 3],
  [0, 1, 0, 3, 2, 3, 2, 1],
  [1, 0, 1, 3, 2, 3, 2, 0],
  [3, 0, 3, 2, 3, 1],
  [3, 0, 3, 2, 3, 1],
  [0, 1, 0, 2, 3, 2, 3, 1],
  [1, 0, 1, 2, 3, 2, 3, 0],
  [2, 0, 2, 1, 2, 3],
  [2, 0, 2, 1, 3, 1, 3, 0],
  [1, 0, 1, 3, 1, 2],
  [0, 1, 0, 2, 0, 3],
  [],
];

function extractSingleSurfaceFromGrid(
  field: Float32Array,
  nx: number,
  ny: number,
  nz: number,
  minX: number,
  minY: number,
  minZ: number,
  dx: number,
  dy: number,
  dz: number,
  isoValue: number,
  phaseSign: 1 | -1
): ColoredIsosurfaceMeshData {
  const vx = nx + 1;
  const vy = ny + 1;
  const vz = nz + 1;
  const slice = vx * vy;
  const totalVerts = vx * vy * vz;

  // Central-difference gradients
  const gradX = new Float32Array(totalVerts);
  const gradY = new Float32Array(totalVerts);
  const gradZ = new Float32Array(totalVerts);
  const inv2Dx = 0.5 / dx;
  const inv2Dy = 0.5 / dy;
  const inv2Dz = 0.5 / dz;

  for (let iz = 0; iz < vz; iz++) {
    const izP = Math.max(0, iz - 1) * slice;
    const izN = Math.min(vz - 1, iz + 1) * slice;
    const zO = iz * slice;
    for (let iy = 0; iy < vy; iy++) {
      const iyP = Math.max(0, iy - 1) * vx;
      const iyN = Math.min(vy - 1, iy + 1) * vx;
      const yO = zO + iy * vx;
      for (let ix = 0; ix < vx; ix++) {
        const idx = yO + ix;
        gradX[idx] = (field[yO + Math.min(vx - 1, ix + 1)] - field[yO + Math.max(0, ix - 1)]) * inv2Dx;
        gradY[idx] = (field[zO + iyN + ix] - field[zO + iyP + ix]) * inv2Dy;
        gradZ[idx] = (field[izN + iy * vx + ix] - field[izP + iy * vx + ix]) * inv2Dz;
      }
    }
  }

  const positions: number[] = [];
  const normals: number[] = [];
  const indices: number[] = [];
  const edgeMap = new Map<number, number>();

  const cOff = [0, 1, 1 + vx, vx, slice, slice + 1, slice + 1 + vx, slice + vx];
  const cX = [0, dx, dx, 0, 0, dx, dx, 0];
  const cY = [0, 0, dy, dy, 0, 0, dy, dy];
  const cZ = [0, 0, 0, 0, dz, dz, dz, dz];

  const getVert = (
    iA: number,
    iB: number,
    ax: number,
    ay: number,
    az: number,
    bx: number,
    by: number,
    bz: number,
    vA: number,
    vB: number
  ) => {
    const lo = iA < iB ? iA : iB;
    const hi = iA < iB ? iB : iA;
    const key = lo * 131072 + hi;
    const found = edgeMap.get(key);
    if (found !== undefined) return found;

    const denom = vB - vA;
    const mu = Math.abs(denom) > 1e-8 ? Math.max(0, Math.min(1, (isoValue - vA) / denom)) : 0.5;
    const px = ax + mu * (bx - ax);
    const py = ay + mu * (by - ay);
    const pz = az + mu * (bz - az);

    let gx = -(gradX[iA] + mu * (gradX[iB] - gradX[iA])) * phaseSign;
    let gy = -(gradY[iA] + mu * (gradY[iB] - gradY[iA])) * phaseSign;
    let gz = -(gradZ[iA] + mu * (gradZ[iB] - gradZ[iA])) * phaseSign;
    const gl = Math.hypot(gx, gy, gz) || 1;
    gx /= gl;
    gy /= gl;
    gz /= gl;

    const vIdx = positions.length / 3;
    positions.push(px, py, pz);
    normals.push(gx, gy, gz);
    edgeMap.set(key, vIdx);
    return vIdx;
  };

  const pushTri = (i0: number, i1: number, i2: number) => {
    const p0x = positions[i0 * 3], p0y = positions[i0 * 3 + 1], p0z = positions[i0 * 3 + 2];
    const ux = positions[i1 * 3] - p0x, uy = positions[i1 * 3 + 1] - p0y, uz = positions[i1 * 3 + 2] - p0z;
    const vx2 = positions[i2 * 3] - p0x, vy2 = positions[i2 * 3 + 1] - p0y, vz2 = positions[i2 * 3 + 2] - p0z;
    const fnx = uy * vz2 - uz * vy2;
    const fny = uz * vx2 - ux * vz2;
    const fnz = ux * vy2 - uy * vx2;
    const gnx = normals[i0 * 3] + normals[i1 * 3] + normals[i2 * 3];
    const gny = normals[i0 * 3 + 1] + normals[i1 * 3 + 1] + normals[i2 * 3 + 1];
    const gnz = normals[i0 * 3 + 2] + normals[i1 * 3 + 2] + normals[i2 * 3 + 2];
    if (fnx * gnx + fny * gny + fnz * gnz >= 0) indices.push(i0, i1, i2);
    else indices.push(i0, i2, i1);
  };

  const cv = new Float32Array(8);
  const ci = new Int32Array(8);
  const cx = new Float32Array(8);
  const cy = new Float32Array(8);
  const cz = new Float32Array(8);

  for (let iz = 0; iz < nz; iz++) {
    const z0 = minZ + iz * dz;
    const zb = iz * slice;
    for (let iy = 0; iy < ny; iy++) {
      const y0 = minY + iy * dy;
      const yb = zb + iy * vx;
      for (let ix = 0; ix < nx; ix++) {
        const base = yb + ix;
        const x0 = minX + ix * dx;
        let above = false;
        let below = false;
        for (let c = 0; c < 8; c++) {
          const gi = base + cOff[c];
          const v = field[gi] * phaseSign;
          cv[c] = v;
          ci[c] = gi;
          if (v > isoValue) above = true;
          else below = true;
        }
        if (!above || !below) continue;
        for (let c = 0; c < 8; c++) {
          cx[c] = x0 + cX[c];
          cy[c] = y0 + cY[c];
          cz[c] = z0 + cZ[c];
        }
        for (let t = 0; t < 6; t++) {
          const tet = KUHN_TETS[t];
          let mask = 0;
          if (cv[tet[0]] > isoValue) mask |= 1;
          if (cv[tet[1]] > isoValue) mask |= 2;
          if (cv[tet[2]] > isoValue) mask |= 4;
          if (cv[tet[3]] > isoValue) mask |= 8;
          if (mask === 0 || mask === 15) continue;
          const ed = TET_EDGES[mask];
          if (ed.length === 6) {
            const i0 = getVert(ci[tet[ed[0]]], ci[tet[ed[1]]], cx[tet[ed[0]]], cy[tet[ed[0]]], cz[tet[ed[0]]], cx[tet[ed[1]]], cy[tet[ed[1]]], cz[tet[ed[1]]], cv[tet[ed[0]]], cv[tet[ed[1]]]);
            const i1 = getVert(ci[tet[ed[2]]], ci[tet[ed[3]]], cx[tet[ed[2]]], cy[tet[ed[2]]], cz[tet[ed[2]]], cx[tet[ed[3]]], cy[tet[ed[3]]], cz[tet[ed[3]]], cv[tet[ed[2]]], cv[tet[ed[3]]]);
            const i2 = getVert(ci[tet[ed[4]]], ci[tet[ed[5]]], cx[tet[ed[4]]], cy[tet[ed[4]]], cz[tet[ed[4]]], cx[tet[ed[5]]], cy[tet[ed[5]]], cz[tet[ed[5]]], cv[tet[ed[4]]], cv[tet[ed[5]]]);
            pushTri(i0, i1, i2);
          } else if (ed.length === 8) {
            const i0 = getVert(ci[tet[ed[0]]], ci[tet[ed[1]]], cx[tet[ed[0]]], cy[tet[ed[0]]], cz[tet[ed[0]]], cx[tet[ed[1]]], cy[tet[ed[1]]], cz[tet[ed[1]]], cv[tet[ed[0]]], cv[tet[ed[1]]]);
            const i1 = getVert(ci[tet[ed[2]]], ci[tet[ed[3]]], cx[tet[ed[2]]], cy[tet[ed[2]]], cz[tet[ed[2]]], cx[tet[ed[3]]], cy[tet[ed[3]]], cz[tet[ed[3]]], cv[tet[ed[2]]], cv[tet[ed[3]]]);
            const i2 = getVert(ci[tet[ed[4]]], ci[tet[ed[5]]], cx[tet[ed[4]]], cy[tet[ed[4]]], cz[tet[ed[4]]], cx[tet[ed[5]]], cy[tet[ed[5]]], cz[tet[ed[5]]], cv[tet[ed[4]]], cv[tet[ed[5]]]);
            const i3 = getVert(ci[tet[ed[6]]], ci[tet[ed[7]]], cx[tet[ed[6]]], cy[tet[ed[6]]], cz[tet[ed[6]]], cx[tet[ed[7]]], cy[tet[ed[7]]], cz[tet[ed[7]]], cv[tet[ed[6]]], cv[tet[ed[7]]]);
            pushTri(i0, i1, i2);
            pushTri(i0, i2, i3);
          }
        }
      }
    }
  }

  return {
    positions: new Float32Array(positions),
    normals: new Float32Array(normals),
    indices: new Uint32Array(indices),
  };
}

/**
 * Converts an sRGB channel in [0, 1] into Linear-sRGB working space for Three.js vertex colors.
 * Without this conversion, Three.js applies an extra gamma lift on vertex attributes, washing
 * saturated red/blue/green ramps into pale pastels.
 */
function srgbChannelToLinear(c: number): number {
  const v = Math.max(0, Math.min(1, c));
  return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4);
}

/**
 * Maps a normalized value t ∈ [-1, 1] onto a high-contrast computational chemistry MEP
 * spectral colormap (sRGB [0..1] output):
 *   -1.00 (nucleophilic, negative V) -> Deep Crimson (#D6080F)
 *   -0.50 -> Vivid Orange (#F77300)
 *   -0.18 -> Rich Gold-Yellow (#F2CC0D)
 *    0.00 (neutral) -> Crisp Emerald (#0FAD47)
 *   +0.18 -> Vivid Turquoise (#00B8C7)
 *   +0.50 -> Rich Azure (#0070F0)
 *   +1.00 (electrophilic, positive V) -> Deep Cobalt (#0F24DB)
 */
export function mepColorRamp(tNorm: number): [number, number, number] {
  const t = Math.max(-1, Math.min(1, tNorm));
  const stops: [number, [number, number, number]][] = [
    [-1.0, [0.84, 0.03, 0.06]],
    [-0.5, [0.97, 0.45, 0.0]],
    [-0.18, [0.95, 0.8, 0.05]],
    [0.0, [0.06, 0.68, 0.28]],
    [0.18, [0.0, 0.72, 0.78]],
    [0.5, [0.0, 0.44, 0.94]],
    [1.0, [0.06, 0.14, 0.86]],
  ];
  for (let i = 0; i < stops.length - 1; i++) {
    const [t0, c0] = stops[i];
    const [t1, c1] = stops[i + 1];
    if (t >= t0 && t <= t1) {
      const f = (t - t0) / (t1 - t0);
      return [
        c0[0] + f * (c1[0] - c0[0]),
        c0[1] + f * (c1[1] - c0[1]),
        c0[2] + f * (c1[2] - c0[2]),
      ];
    }
  }
  return [0.06, 0.68, 0.28];
}

/**
 * Builds 3D isosurfaces for Total Electron Density (ρ), MEP-on-Density (V_ESP on ρ),
 * Spin Density (ρ_α - ρ_β), or an imported Gaussian .cube file.
 */
export function buildVolumetricSurfaceV3(params: {
  fieldMode: Exclude<VolumetricFieldMode, 'mo'>;
  atoms: AtomNode[];
  eht: LowdinEHTResult;
  lobeScale: number;
  importedCube?: GaussianCubeData | null;
}): VolumetricSurfaceResult {
  const { fieldMode, atoms, eht, lobeScale, importedCube } = params;
  const empty: ColoredIsosurfaceMeshData = {
    positions: new Float32Array(0),
    normals: new Float32Array(0),
    indices: new Uint32Array(0),
  };

  if (fieldMode === 'cube') {
    if (!importedCube) return { positive: empty, negative: empty, hasData: false };
    const nx = Math.min(36, importedCube.nx - 1);
    const ny = Math.min(36, importedCube.ny - 1);
    const nz = Math.min(36, importedCube.nz - 1);
    const spanX = (importedCube.nx - 1) * (importedCube.vxAngstrom[0] || 0.2);
    const spanY = (importedCube.ny - 1) * (importedCube.vyAngstrom[1] || 0.2);
    const spanZ = (importedCube.nz - 1) * (importedCube.vzAngstrom[2] || 0.2);
    const dx = spanX / nx;
    const dy = spanY / ny;
    const dz = spanZ / nz;
    const vx = nx + 1, vy = ny + 1, vz = nz + 1;
    const field = new Float32Array(vx * vy * vz);
    const slice = vx * vy;
    let absMax = 0;

    for (let iz = 0; iz < vz; iz++) {
      const z = importedCube.originAngstrom[2] + iz * dz;
      const wz = iz === 0 || iz === nz ? 0 : 1;
      for (let iy = 0; iy < vy; iy++) {
        const y = importedCube.originAngstrom[1] + iy * dy;
        const wy = iy === 0 || iy === ny ? 0 : wz;
        for (let ix = 0; ix < vx; ix++) {
          if (wy === 0 || ix === 0 || ix === nx) continue;
          const x = importedCube.originAngstrom[0] + ix * dx;
          const val = sampleGaussianCubeAtPoint(importedCube, x, y, z);
          field[iz * slice + iy * vx + ix] = val;
          if (Math.abs(val) > absMax) absMax = Math.abs(val);
        }
      }
    }
    const iso = Math.max(1e-4, absMax * 0.14) / Math.max(0.45, Math.min(1.65, lobeScale));
    const positive = extractSingleSurfaceFromGrid(field, nx, ny, nz, importedCube.originAngstrom[0], importedCube.originAngstrom[1], importedCube.originAngstrom[2], dx, dy, dz, iso, 1);
    const hasNeg = importedCube.minVal < -absMax * 0.05;
    const negative = hasNeg
      ? extractSingleSurfaceFromGrid(field, nx, ny, nz, importedCube.originAngstrom[0], importedCube.originAngstrom[1], importedCube.originAngstrom[2], dx, dy, dz, iso, -1)
      : empty;
    return {
      positive,
      negative,
      hasData: positive.positions.length > 0 || negative.positions.length > 0,
    };
  }

  if (atoms.length === 0) return { positive: empty, negative: empty, hasData: false };
  if (fieldMode === 'spin' && eht.unpairedElectrons === 0) {
    return { positive: empty, negative: empty, hasData: false };
  }

  let minX = Infinity, maxX = -Infinity;
  let minY = Infinity, maxY = -Infinity;
  let minZ = Infinity, maxZ = -Infinity;
  for (const a of atoms) {
    if (a.x3d < minX) minX = a.x3d;
    if (a.x3d > maxX) maxX = a.x3d;
    if (a.y3d < minY) minY = a.y3d;
    if (a.y3d > maxY) maxY = a.y3d;
    if (a.z3d < minZ) minZ = a.z3d;
    if (a.z3d > maxZ) maxZ = a.z3d;
  }
  const pad = fieldMode === 'mep' ? 2.9 : 2.6;
  minX -= pad; maxX += pad;
  minY -= pad; maxY += pad;
  minZ -= pad; maxZ += pad;
  const spanX = Math.max(4.4, maxX - minX);
  const spanY = Math.max(4.4, maxY - minY);
  const spanZ = Math.max(4.4, maxZ - minZ);
  const maxSpan = Math.max(spanX, spanY, spanZ);

  const maxCells = atoms.length > 14 ? 26 : fieldMode === 'mep' ? 32 : 30;
  const nx = Math.max(18, Math.round((spanX / maxSpan) * maxCells));
  const ny = Math.max(18, Math.round((spanY / maxSpan) * maxCells));
  const nz = Math.max(18, Math.round((spanZ / maxSpan) * maxCells));
  const dx = spanX / nx;
  const dy = spanY / ny;
  const dz = spanZ / nz;
  const vx = nx + 1, vy = ny + 1, vz = nz + 1;
  const slice = vx * vy;
  const field = new Float32Array(vx * vy * vz);

  let maxAbs = 0;
  for (let iz = 0; iz < vz; iz++) {
    const z = minZ + iz * dz;
    const wz = iz === 0 || iz === nz ? 0 : 1;
    for (let iy = 0; iy < vy; iy++) {
      const y = minY + iy * dy;
      const wy = iy === 0 || iy === ny ? 0 : wz;
      for (let ix = 0; ix < vx; ix++) {
        if (wy === 0 || ix === 0 || ix === nx) continue;
        const x = minX + ix * dx;
        const res = evalDensityAndSpinAtPoint(eht, atoms, x, y, z);
        const val = fieldMode === 'spin' ? res.spinDensity : res.density;
        field[iz * slice + iy * vx + ix] = val;
        if (Math.abs(val) > maxAbs) maxAbs = Math.abs(val);
      }
    }
  }

  if (maxAbs < 1e-7) return { positive: empty, negative: empty, hasData: false };

  // Calibrated chemical isovalues:
  //   - MEP: van der Waals outer density envelope (~0.003-0.005 e/a0^3, standard Gaussian/Spartan convention)
  //   - Density: 0.02-0.04 e/a0^3 covalent envelope
  //   - Spin: 0.002+ e/a0^3 spin excess lobe
  const baseIso =
    fieldMode === 'spin'
      ? Math.max(0.002, maxAbs * 0.08)
      : fieldMode === 'mep'
        ? Math.max(0.0028, Math.min(0.0052, maxAbs * 0.008))
        : Math.max(0.012, Math.min(0.045, maxAbs * 0.055));
  const iso = baseIso / Math.max(0.45, Math.min(1.65, lobeScale));

  const positive = extractSingleSurfaceFromGrid(field, nx, ny, nz, minX, minY, minZ, dx, dy, dz, iso, 1);

  if (fieldMode === 'mep' && positive.positions.length > 0) {
    const vCount = positive.positions.length / 3;
    const potentials = new Float32Array(vCount);
    const absPotentials = new Float32Array(vCount);
    for (let i = 0; i < vCount; i++) {
      const v = evalElectrostaticPotentialAtPoint(
        eht,
        atoms,
        positive.positions[i * 3],
        positive.positions[i * 3 + 1],
        positive.positions[i * 3 + 2]
      );
      potentials[i] = v;
      absPotentials[i] = Math.abs(v);
    }
    // Use 94th-percentile symmetric bound (floored at 0.015 a.u.) so both polar
    // heterocycles/water and quadrupolar aromatics/N2 render with crisp, faithful contrast.
    absPotentials.sort();
    const p94 = absPotentials[Math.min(vCount - 1, Math.floor(vCount * 0.94))] ?? 0;
    const vScale = Math.max(0.015, p94);

    const colors = new Float32Array(vCount * 3);
    for (let i = 0; i < vCount; i++) {
      const raw = Math.max(-1, Math.min(1, potentials[i] / vScale));
      const shaped = Math.sign(raw) * Math.pow(Math.abs(raw), 0.85);
      const [rSrgb, gSrgb, bSrgb] = mepColorRamp(shaped);
      colors[i * 3] = srgbChannelToLinear(rSrgb);
      colors[i * 3 + 1] = srgbChannelToLinear(gSrgb);
      colors[i * 3 + 2] = srgbChannelToLinear(bSrgb);
    }
    positive.colors = colors;
    return {
      positive,
      negative: empty,
      mepRange: { min: -vScale, max: vScale },
      hasData: true,
    };
  }

  const negative =
    fieldMode === 'spin'
      ? extractSingleSurfaceFromGrid(field, nx, ny, nz, minX, minY, minZ, dx, dy, dz, iso, -1)
      : empty;

  return {
    positive,
    negative,
    hasData: positive.positions.length > 0 || negative.positions.length > 0,
  };
}

/**
 * Evaluates a uniform 3D grid (x slowest, y middle, z fastest) for Gaussian .cube export.
 */
export function evaluateUniformGridForCubeExport(params: {
  fieldMode: VolumetricFieldMode;
  atoms: AtomNode[];
  selectedMo: MolecularOrbital;
  eht: LowdinEHTResult;
  importedCube?: GaussianCubeData | null;
  resolution?: number;
}): {
  originAngstrom: [number, number, number];
  nx: number;
  ny: number;
  nz: number;
  dxAngstrom: number;
  dyAngstrom: number;
  dzAngstrom: number;
  values: Float32Array;
} {
  const { fieldMode, atoms, selectedMo, eht, importedCube, resolution = 24 } = params;
  if (fieldMode === 'cube' && importedCube) {
    return {
      originAngstrom: importedCube.originAngstrom,
      nx: importedCube.nx,
      ny: importedCube.ny,
      nz: importedCube.nz,
      dxAngstrom: importedCube.vxAngstrom[0] || 0.2,
      dyAngstrom: importedCube.vyAngstrom[1] || 0.2,
      dzAngstrom: importedCube.vzAngstrom[2] || 0.2,
      values: importedCube.values,
    };
  }

  let minX = -2, maxX = 2, minY = -2, maxY = 2, minZ = -2, maxZ = 2;
  if (atoms.length > 0) {
    minX = Math.min(...atoms.map((a) => a.x3d)) - 2.6;
    maxX = Math.max(...atoms.map((a) => a.x3d)) + 2.6;
    minY = Math.min(...atoms.map((a) => a.y3d)) - 2.6;
    maxY = Math.max(...atoms.map((a) => a.y3d)) + 2.6;
    minZ = Math.min(...atoms.map((a) => a.z3d)) - 2.6;
    maxZ = Math.max(...atoms.map((a) => a.z3d)) + 2.6;
  }
  const nx = resolution;
  const ny = resolution;
  const nz = resolution;
  const dxAngstrom = (maxX - minX) / (nx - 1);
  const dyAngstrom = (maxY - minY) / (ny - 1);
  const dzAngstrom = (maxZ - minZ) / (nz - 1);
  const values = new Float32Array(nx * ny * nz);

  for (let ix = 0; ix < nx; ix++) {
    const x = minX + ix * dxAngstrom;
    for (let iy = 0; iy < ny; iy++) {
      const y = minY + iy * dyAngstrom;
      for (let iz = 0; iz < nz; iz++) {
        const z = minZ + iz * dzAngstrom;
        let v = 0;
        if (fieldMode === 'mo') {
          v = evalRealSTOWavefunction(selectedMo.lcao, atoms, x, y, z);
        } else if (fieldMode === 'mep') {
          v = evalElectrostaticPotentialAtPoint(eht, atoms, x, y, z);
        } else {
          const ds = evalDensityAndSpinAtPoint(eht, atoms, x, y, z);
          v = fieldMode === 'spin' ? ds.spinDensity : ds.density;
        }
        values[(ix * ny + iy) * nz + iz] = v;
      }
    }
  }

  return {
    originAngstrom: [minX, minY, minZ],
    nx,
    ny,
    nz,
    dxAngstrom,
    dyAngstrom,
    dzAngstrom,
    values,
  };
}

/**
 * Computes a crisp 2D contour & diverging-phase slice plane across xy, xz, or yz at `offset` Angstroms.
 */
export function compute2DSlicePlane(params: {
  fieldMode: VolumetricFieldMode;
  atoms: AtomNode[];
  selectedMo: MolecularOrbital;
  eht: LowdinEHTResult;
  orientation: SliceOrientation;
  offset: number;
  showIsolines: boolean;
  colorblindPalette?: boolean;
  importedCube?: GaussianCubeData | null;
  resolution?: number;
}): Slice2DResult {
  const {
    fieldMode,
    atoms,
    selectedMo,
    eht,
    orientation,
    offset,
    showIsolines,
    colorblindPalette = false,
    importedCube,
    resolution = 84,
  } = params;

  const rMol = Math.max(1.8, ...atoms.map((a) => Math.hypot(a.x3d, a.y3d, a.z3d)));
  const halfSpan = Math.min(8.5, rMol + 2.2);
  const spanU = halfSpan * 2;
  const spanV = halfSpan * 2;

  const width = resolution;
  const height = resolution;
  const values = new Float32Array(width * height);
  let minVal = Infinity;
  let maxVal = -Infinity;

  for (let iy = 0; iy < height; iy++) {
    const vCoord = -halfSpan + (iy / (height - 1)) * spanV;
    for (let ix = 0; ix < width; ix++) {
      const uCoord = -halfSpan + (ix / (width - 1)) * spanU;
      let x = 0, y = 0, z = 0;
      if (orientation === 'xy') {
        x = uCoord; y = vCoord; z = offset;
      } else if (orientation === 'xz') {
        x = uCoord; y = offset; z = vCoord;
      } else {
        x = offset; y = vCoord; z = uCoord;
      }

      let val = 0;
      if (fieldMode === 'mo') {
        val = evalRealSTOWavefunction(selectedMo.lcao, atoms, x, y, z);
      } else if (fieldMode === 'mep') {
        val = evalElectrostaticPotentialAtPoint(eht, atoms, x, y, z);
      } else if (fieldMode === 'cube' && importedCube) {
        val = sampleGaussianCubeAtPoint(importedCube, x, y, z);
      } else {
        const ds = evalDensityAndSpinAtPoint(eht, atoms, x, y, z);
        val = fieldMode === 'spin' ? ds.spinDensity : ds.density;
      }
      values[iy * width + ix] = val;
      if (val < minVal) minVal = val;
      if (val > maxVal) maxVal = val;
    }
  }

  const absMax = Math.max(1e-5, Math.abs(minVal), Math.abs(maxVal));
  let mepSliceScale = absMax;
  if (fieldMode === 'mep') {
    const absArr = new Float32Array(values.length);
    for (let i = 0; i < values.length; i++) absArr[i] = Math.abs(values[i]);
    absArr.sort();
    const p92 = absArr[Math.min(absArr.length - 1, Math.floor(absArr.length * 0.92))] ?? 0;
    mepSliceScale = Math.max(0.02, p92);
  }
  const rgba = new Uint8ClampedArray(width * height * 4);

  // Contour levels (logarithmic-like spacing for smooth orbital contours)
  const levels = [0.06, 0.14, 0.28, 0.50, 0.78];

  for (let iy = 0; iy < height; iy++) {
    for (let ix = 0; ix < width; ix++) {
      const val = values[iy * width + ix];
      const norm = val / (fieldMode === 'mep' ? mepSliceScale : absMax);
      const mag = Math.min(1, Math.pow(Math.abs(norm), 0.55));

      let r = 248, g = 250, b = 252, a = 195;

      if (fieldMode === 'mep') {
        const clamped = Math.max(-1, Math.min(1, norm));
        const shaped = Math.sign(clamped) * Math.pow(Math.abs(clamped), 0.85);
        const [mr, mg, mb] = mepColorRamp(shaped);
        r = Math.round(mr * 255);
        g = Math.round(mg * 255);
        b = Math.round(mb * 255);
        a = 210;
      } else if (fieldMode === 'density') {
        // Viridis-inspired sequential density map (transparent near 0 -> teal -> bright amber)
        const t = Math.min(1, Math.pow(Math.max(0, norm), 0.45));
        r = Math.round(15 + t * 235);
        g = Math.round(23 + Math.sin(t * Math.PI * 0.85) * 190);
        b = Math.round(60 + (1 - t) * 175);
        a = Math.round(40 + t * 195);
      } else {
        // Signed diverging phase colormap (Cobalt/Vermilion or Teal/Amber)
        const posRGB = colorblindPalette ? [13, 148, 136] : [37, 99, 235];
        const negRGB = colorblindPalette ? [217, 119, 6] : [220, 38, 38];
        const target = val >= 0 ? posRGB : negRGB;
        r = Math.round(248 * (1 - mag) + target[0] * mag);
        g = Math.round(250 * (1 - mag) + target[1] * mag);
        b = Math.round(252 * (1 - mag) + target[2] * mag);
        a = Math.round(55 + mag * 180);
      }

      // Crisp isoline contour stroke detection (zero nodal line + signed amplitude levels)
      if (showIsolines && ix < width - 1 && iy < height - 1) {
        const vRight = values[iy * width + (ix + 1)] / absMax;
        const vDown = values[(iy + 1) * width + ix] / absMax;
        const vSelf = norm;
        const isNodeCross =
          fieldMode !== 'density' &&
          ((vSelf >= 0 && (vRight < 0 || vDown < 0)) || (vSelf < 0 && (vRight >= 0 || vDown >= 0)));
        if (isNodeCross) {
          r = 15; g = 23; b = 42; a = 235;
        } else {
          for (const lv of levels) {
            const cPos = (vSelf - lv) * (vRight - lv) <= 0 || (vSelf - lv) * (vDown - lv) <= 0;
            const cNeg =
              fieldMode !== 'density' &&
              ((vSelf + lv) * (vRight + lv) <= 0 || (vSelf + lv) * (vDown + lv) <= 0);
            if (cPos || cNeg) {
              r = Math.round(r * 0.45);
              g = Math.round(g * 0.45);
              b = Math.round(b * 0.45);
              a = 240;
              break;
            }
          }
        }
      }

      const pIdx = (iy * width + ix) * 4;
      rgba[pIdx] = r;
      rgba[pIdx + 1] = g;
      rgba[pIdx + 2] = b;
      rgba[pIdx + 3] = a;
    }
  }

  const center: [number, number, number] =
    orientation === 'xy' ? [0, 0, offset] : orientation === 'xz' ? [0, offset, 0] : [offset, 0, 0];

  return {
    width,
    height,
    values,
    rgba,
    minVal: Number.isFinite(minVal) ? minVal : 0,
    maxVal: Number.isFinite(maxVal) ? maxVal : 0,
    absMax,
    spanU,
    spanV,
    center,
    orientation,
    offset,
  };
}
