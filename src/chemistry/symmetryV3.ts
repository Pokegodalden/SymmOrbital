import {
  AtomNode,
  ELEMENT_PROPERTIES,
  MolecularOrbital,
  OrbitalBasisType,
} from './moleculesData';
import {
  POINT_GROUPS,
  PointGroupData,
  SymmetryElement3D,
  decomposeReducibleRepresentation,
} from './symmetryData';
import {
  LowdinEHTResult,
  STOBasisFunction,
  symmetricJacobiEigen,
} from './quantumEngineV3';

export type Mat3 = [
  [number, number, number],
  [number, number, number],
  [number, number, number],
];
export type Vec3 = [number, number, number];

export interface SymmetryDetectionV3 {
  pointGroupId: string;
  pointGroup: PointGroupData;
  centroid: Vec3;
  principalAxes: Mat3; // Columns are orthonormal principal axes [eX, eY, eZ] in molecular frame
  principalMoments: Vec3;
  alignedElements3D: SymmetryElement3D[];
  hasInversion: boolean;
  principalFold: number;
  rmsdTolerance: number;
}

export interface MOIrrepProjection {
  moIndex: number;
  irrepId: string;
  irrepSymbol: string;
  confidence: number;
  computedCharacters: number[];
}

export interface ElectronicTransitionRule {
  fromMoId: string;
  fromLabel: string;
  fromIrrepId: string;
  fromIrrepSymbol: string;
  toMoId: string;
  toLabel: string;
  toIrrepId: string;
  toIrrepSymbol: string;
  deltaEv: number;
  wavelengthNm: number;
  wavenumberCm1: number;
  directProductFormula: string;
  allowedPolarizations: ('x' | 'y' | 'z')[];
  isElectricDipoleAllowed: boolean;
  laporteStatus: 'allowed' | 'forbidden' | 'not-applicable';
  transitionDipoleStrength: number; // Approximate |<i|r|f>| in Angstrom
  rationale: string;
}

export interface VibrationalActivitySummary {
  isCentrosymmetric: boolean;
  mutualExclusionApplies: boolean;
  irActiveIrreps: string[];
  ramanActiveIrreps: string[];
  coincidentIrreps: string[];
  silentIrreps: string[];
}

// ============================================================================
// 3D Linear Algebra & Symmetry Operation Helpers
// ============================================================================

function dot3(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cross3(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function norm3(a: Vec3): number {
  return Math.hypot(a[0], a[1], a[2]);
}

function normalize3(a: Vec3): Vec3 {
  const n = norm3(a);
  if (n < 1e-12) return [0, 0, 1];
  return [a[0] / n, a[1] / n, a[2] / n];
}

function mat3MulVec(M: Mat3, v: Vec3): Vec3 {
  return [
    M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
    M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
    M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
  ];
}

function mat3Mul(A: Mat3, B: Mat3): Mat3 {
  const C: Mat3 = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
    }
  }
  return C;
}

export function rotationMatrixAroundAxis(axis: Vec3, angleRad: number): Mat3 {
  const [ux, uy, uz] = normalize3(axis);
  const c = Math.cos(angleRad);
  const s = Math.sin(angleRad);
  const t = 1 - c;
  return [
    [t * ux * ux + c, t * ux * uy - s * uz, t * ux * uz + s * uy],
    [t * ux * uy + s * uz, t * uy * uy + c, t * uy * uz - s * ux],
    [t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c],
  ];
}

export function reflectionMatrixAcrossPlane(normal: Vec3): Mat3 {
  const [nx, ny, nz] = normalize3(normal);
  return [
    [1 - 2 * nx * nx, -2 * nx * ny, -2 * nx * nz],
    [-2 * nx * ny, 1 - 2 * ny * ny, -2 * ny * nz],
    [-2 * nx * nz, -2 * ny * nz, 1 - 2 * nz * nz],
  ];
}

export function improperRotationMatrix(axis: Vec3, angleRad: number): Mat3 {
  return mat3Mul(reflectionMatrixAcrossPlane(axis), rotationMatrixAroundAxis(axis, angleRad));
}

const INVERSION_MAT: Mat3 = [
  [-1, 0, 0],
  [0, -1, 0],
  [0, 0, -1],
];

const IDENTITY_MAT: Mat3 = [
  [1, 0, 0],
  [0, 1, 0],
  [0, 0, 1],
];

/**
 * Checks whether orthogonal 3x3 matrix R maps every atom to an identical-element atom within `tol` Angstrom.
 */
function checkSymmetryOperation(
  coords: Vec3[],
  elements: string[],
  R: Mat3,
  tol: number,
  requireNonTrivialPermutation = false
): boolean {
  const n = coords.length;
  let movedAnyAtom = false;
  for (let i = 0; i < n; i++) {
    const p = mat3MulVec(R, coords[i]);
    let matchedIdx = -1;
    let bestD = Infinity;
    for (let j = 0; j < n; j++) {
      if (elements[j] !== elements[i]) continue;
      const d = Math.hypot(p[0] - coords[j][0], p[1] - coords[j][1], p[2] - coords[j][2]);
      if (d <= tol && d < bestD) {
        bestD = d;
        matchedIdx = j;
      }
    }
    if (matchedIdx < 0) return false;
    if (matchedIdx !== i) movedAnyAtom = true;
  }
  if (requireNonTrivialPermutation && !movedAnyAtom) return false;
  return true;
}

function dedupUnitVectors(vecs: Vec3[], cosTol = 0.992): Vec3[] {
  const out: Vec3[] = [];
  for (const raw of vecs) {
    const n = norm3(raw);
    if (n < 1e-5) continue;
    const u: Vec3 = [raw[0] / n, raw[1] / n, raw[2] / n];
    if (!out.some((v) => Math.abs(dot3(u, v)) > cosTol)) {
      out.push(u);
    }
  }
  return out;
}

/**
 * Generates candidate symmetry axes and mirror plane normals from principal axes,
 * atom vectors, pair midpoints, pair differences, and triplet normals.
 */
function generateCandidateAxes(coords: Vec3[], elements: string[], principalAxes: Vec3[]): Vec3[] {
  const cands: Vec3[] = [...principalAxes, [1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const n = coords.length;

  for (let i = 0; i < n; i++) {
    if (norm3(coords[i]) > 0.1) cands.push(coords[i]);
  }

  // Limit pairwise sampling for very large macromolecules while staying complete for N <= 40
  const maxN = Math.min(n, 48);
  for (let i = 0; i < maxN; i++) {
    for (let j = i + 1; j < maxN; j++) {
      if (elements[i] !== elements[j]) continue;
      const mid: Vec3 = [
        coords[i][0] + coords[j][0],
        coords[i][1] + coords[j][1],
        coords[i][2] + coords[j][2],
      ];
      if (norm3(mid) > 0.12) cands.push(mid);
      const diff: Vec3 = [
        coords[i][0] - coords[j][0],
        coords[i][1] - coords[j][1],
        coords[i][2] - coords[j][2],
      ];
      if (norm3(diff) > 0.12) cands.push(diff);
      const cr = cross3(coords[i], coords[j]);
      if (norm3(cr) > 0.12) cands.push(cr);
    }
  }

  // Triplet normals for rings & polyhedral faces (small-to-medium systems)
  const triN = Math.min(n, 24);
  for (let i = 0; i < triN; i++) {
    for (let j = i + 1; j < triN; j++) {
      if (elements[i] !== elements[j]) continue;
      for (let k = j + 1; k < triN; k++) {
        if (elements[i] !== elements[k]) continue;
        const a: Vec3 = [
          coords[j][0] - coords[i][0],
          coords[j][1] - coords[i][1],
          coords[j][2] - coords[i][2],
        ];
        const b: Vec3 = [
          coords[k][0] - coords[i][0],
          coords[k][1] - coords[i][1],
          coords[k][2] - coords[i][2],
        ];
        const cr = cross3(a, b);
        if (norm3(cr) > 0.12) cands.push(cr);
      }
    }
  }

  return dedupUnitVectors(cands, 0.9995);
}

/**
 * Exact 27-group 3D point-group detection using mass-weighted principal inertia axes
 * and full 3D operation testing, emitting symmetry elements rotated into the molecule's frame.
 */
export function detectExactPointGroup3D(
  atoms: AtomNode[],
  tolerance = 0.08
): SymmetryDetectionV3 {
  if (atoms.length === 0) {
    const pg = POINT_GROUPS.C1;
    return {
      pointGroupId: 'C1',
      pointGroup: pg,
      centroid: [0, 0, 0],
      principalAxes: IDENTITY_MAT,
      principalMoments: [0, 0, 0],
      alignedElements3D: [],
      hasInversion: false,
      principalFold: 1,
      rmsdTolerance: tolerance,
    };
  }

  // 1. Center of mass
  let totalMass = 0;
  let cmx = 0,
    cmy = 0,
    cmz = 0;
  const masses = atoms.map((a) => {
    const p = ELEMENT_PROPERTIES[a.element];
    return p ? Math.max(1, p.atomicNumber * 2) : 12;
  });
  for (let i = 0; i < atoms.length; i++) {
    totalMass += masses[i];
    cmx += masses[i] * atoms[i].x3d;
    cmy += masses[i] * atoms[i].y3d;
    cmz += masses[i] * atoms[i].z3d;
  }
  const centroid: Vec3 = [cmx / totalMass, cmy / totalMass, cmz / totalMass];
  const coords: Vec3[] = atoms.map((a) => [
    a.x3d - centroid[0],
    a.y3d - centroid[1],
    a.z3d - centroid[2],
  ]);
  const elements = atoms.map((a) => a.element);

  // 2. Mass-weighted inertia tensor
  const I: number[][] = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
  ];
  for (let i = 0; i < coords.length; i++) {
    const [x, y, z] = coords[i];
    const m = masses[i];
    I[0][0] += m * (y * y + z * z);
    I[1][1] += m * (x * x + z * z);
    I[2][2] += m * (x * x + y * y);
    I[0][1] -= m * x * y;
    I[0][2] -= m * x * z;
    I[1][2] -= m * y * z;
  }
  I[1][0] = I[0][1];
  I[2][0] = I[0][2];
  I[2][1] = I[1][2];

  const eig = symmetricJacobiEigen(I);
  const moments: Vec3 = [eig.eigenvalues[0], eig.eigenvalues[1], eig.eigenvalues[2]];
  let e0: Vec3 = normalize3([eig.eigenvectors[0][0], eig.eigenvectors[1][0], eig.eigenvectors[2][0]]);
  let e1: Vec3 = normalize3([eig.eigenvectors[0][1], eig.eigenvectors[1][1], eig.eigenvectors[2][1]]);
  let e2: Vec3 = normalize3(cross3(e0, e1)); // Enforce right-handed frame

  const hasInversion = checkSymmetryOperation(coords, elements, INVERSION_MAT, tolerance);

  // 3. Linear molecule check (Ia ~ 0, Ib ≈ Ic > 0)
  if (atoms.length === 1) {
    return buildResult('Oh', centroid, [e0, e1, e2], moments, true, 4, tolerance);
  }
  const maxPerpDistFromE0 = coords.reduce((mx, p) => {
    const perp = norm3(cross3(p, e0));
    return perp > mx ? perp : mx;
  }, 0);
  const isLinear =
    atoms.length === 2 ||
    (maxPerpDistFromE0 < 0.05 &&
      moments[0] < 2e-3 * Math.max(1, moments[2]) &&
      Math.abs(moments[1] - moments[2]) < 0.15 * Math.max(1, moments[2]));
  if (isLinear) {
    const pgId = hasInversion ? 'Dinfh' : 'Cinfv';
    // Put molecular axis along eZ column
    const linAxis = e0;
    const perp1 = e1;
    const perp2 = normalize3(cross3(linAxis, perp1));
    return buildResult(pgId, centroid, [perp1, perp2, linAxis], moments, hasInversion, 99, tolerance, {
      principalAxis: linAxis,
      perpC2: hasInversion ? perp1 : undefined,
      sigmaH: hasInversion ? linAxis : undefined,
      sigmaVs: [perp1, perp2],
    });
  }

  // 4. Generate candidate axes & mirror planes
  const cands = generateCandidateAxes(coords, elements, [e0, e1, e2]);
  const axesByFold: Record<number, Vec3[]> = { 2: [], 3: [], 4: [], 5: [], 6: [] };
  const mirrors: Vec3[] = [];

  for (const u of cands) {
    if (checkSymmetryOperation(coords, elements, reflectionMatrixAcrossPlane(u), tolerance)) {
      mirrors.push(u);
    }
    for (const fold of [6, 5, 4, 3, 2]) {
      const R = rotationMatrixAroundAxis(u, (2 * Math.PI) / fold);
      // For any non-linear molecule, a genuine C_n axis must permute at least one pair of off-axis atoms
      if (checkSymmetryOperation(coords, elements, R, tolerance, true)) {
        axesByFold[fold].push(u);
      }
    }
  }

  const c5 = dedupUnitVectors(axesByFold[5]);
  const c4 = dedupUnitVectors(axesByFold[4]);
  const c3 = dedupUnitVectors(axesByFold[3]);
  const c6 = dedupUnitVectors(axesByFold[6]);
  const c2 = dedupUnitVectors(axesByFold[2]);
  const sigmaList = dedupUnitVectors(mirrors);

  // 5. Cubic / Polyhedral groups (multiple C3/C4/C5 axes)
  if (c5.length >= 6 && hasInversion) {
    return buildResult('Ih', centroid, [e0, e1, e2], moments, true, 5, tolerance, {
      principalAxis: c5[0],
      extraAxes: [c3[0]].filter(Boolean) as Vec3[],
    });
  }
  if (c4.length >= 3 && c3.length >= 4 && hasInversion) {
    return buildResult('Oh', centroid, [c4[0], c4[1], c4[2]], moments, true, 4, tolerance, {
      principalAxis: c4[2] || e2,
      extraAxes: [c4[0], c4[1], c3[0]].filter(Boolean) as Vec3[],
      sigmaH: c4[2],
      sigmaVs: [c4[0], c4[1]].filter(Boolean) as Vec3[],
    });
  }
  if (c3.length >= 4) {
    if (hasInversion) {
      return buildResult('Th', centroid, [c2[0] || e0, c2[1] || e1, c2[2] || e2], moments, true, 3, tolerance, {
        principalAxis: c2[2] || e2,
        extraAxes: [c3[0]],
        sigmaH: c2[2] || e2,
      });
    }
    if (sigmaList.length >= 6) {
      return buildResult('Td', centroid, [c2[0] || e0, c2[1] || e1, c2[2] || e2], moments, false, 3, tolerance, {
        principalAxis: c3[0],
        extraAxes: [c2[0] || e2],
        sigmaVs: sigmaList.slice(0, 2),
      });
    }
  }

  // 6. Determine highest-order principal rotation axis
  let principalFold = 1;
  let principalAxis: Vec3 | null = null;
  for (const f of [6, 5, 4, 3, 2]) {
    if (axesByFold[f].length > 0) {
      principalFold = f;
      // For D2/D2h/C2v with multiple C2 axes, prefer the C2 axis with most atoms on it or unique inertia axis
      if (f === 2 && c2.length > 1) {
        let bestAxis = c2[0];
        let bestScore = -1;
        for (const ax of c2) {
          let onAxis = 0;
          for (const p of coords) {
            if (norm3(cross3(p, ax)) < tolerance) onAxis++;
          }
          // Slightly prefer axis aligned with z or y for readable standard orientation
          const score = onAxis * 10 + Math.abs(ax[2]) * 0.5 + Math.abs(ax[1]) * 0.25;
          if (score > bestScore) {
            bestScore = score;
            bestAxis = ax;
          }
        }
        principalAxis = bestAxis;
      } else {
        principalAxis = axesByFold[f][0];
      }
      break;
    }
  }

  // 7. Classify non-axial (C1, Cs, Ci)
  if (!principalAxis || principalFold === 1) {
    if (hasInversion) {
      return buildResult('Ci', centroid, [e0, e1, e2], moments, true, 1, tolerance);
    }
    if (sigmaList.length > 0) {
      const nPlane = sigmaList[0];
      const u1 = normalize3(cross3(nPlane, Math.abs(nPlane[0]) < 0.8 ? [1, 0, 0] : [0, 1, 0]));
      const u2 = normalize3(cross3(nPlane, u1));
      return buildResult('Cs', centroid, [u1, u2, nPlane], moments, false, 1, tolerance, {
        sigmaH: nPlane,
      });
    }
    return buildResult('C1', centroid, [e0, e1, e2], moments, false, 1, tolerance);
  }

  // 8. Axial groups (Cn, Cnv, Cnh, Sn, Dn, Dnh, Dnd)
  // Explicitly probe exact in-plane perpendicular C2' axes and σ_v/σ_d mirror normals around principalAxis
  // so large symmetric molecules (e.g. Decamethylferrocene D5d, Coronene D6h) never lose a true C2'/σ_d
  // axis to off-plane candidate contamination.
  const inPlaneCandsRaw: Vec3[] = [];
  const perpSeedA = normalize3(
    cross3(principalAxis, Math.abs(principalAxis[0]) < 0.8 ? [1, 0, 0] : [0, 1, 0])
  );
  const perpSeedB = normalize3(cross3(principalAxis, perpSeedA));
  inPlaneCandsRaw.push(perpSeedA, perpSeedB);
  const projPerp: Vec3[] = [];
  for (let i = 0; i < coords.length; i++) {
    const p = coords[i];
    const dAxis = dot3(p, principalAxis);
    const v: Vec3 = [
      p[0] - dAxis * principalAxis[0],
      p[1] - dAxis * principalAxis[1],
      p[2] - dAxis * principalAxis[2],
    ];
    const vLen = norm3(v);
    if (vLen > 0.15) {
      const u: Vec3 = [v[0] / vLen, v[1] / vLen, v[2] / vLen];
      projPerp.push(u);
      inPlaneCandsRaw.push(u, cross3(principalAxis, u));
    }
  }
  const uniqueProj = dedupUnitVectors(projPerp, 0.9995);
  for (let i = 0; i < Math.min(uniqueProj.length, 32); i++) {
    for (let j = i + 1; j < Math.min(uniqueProj.length, 32); j++) {
      const sVec: Vec3 = [
        uniqueProj[i][0] + uniqueProj[j][0],
        uniqueProj[i][1] + uniqueProj[j][1],
        uniqueProj[i][2] + uniqueProj[j][2],
      ];
      if (norm3(sVec) > 0.1) inPlaneCandsRaw.push(sVec, cross3(principalAxis, sVec));
      const dVec: Vec3 = [
        uniqueProj[i][0] - uniqueProj[j][0],
        uniqueProj[i][1] - uniqueProj[j][1],
        uniqueProj[i][2] - uniqueProj[j][2],
      ];
      if (norm3(dVec) > 0.1) inPlaneCandsRaw.push(dVec, cross3(principalAxis, dVec));
    }
  }
  const inPlaneCands = dedupUnitVectors(inPlaneCandsRaw, 0.9995);
  const exactPerpC2Raw: Vec3[] = [...c2.filter((ax) => Math.abs(dot3(ax, principalAxis!)) < 0.12)];
  const exactSigmaVRaw: Vec3[] = [
    ...sigmaList.filter((m) => Math.abs(dot3(m, principalAxis!)) < 0.12),
  ];
  for (const u of inPlaneCands) {
    if (checkSymmetryOperation(coords, elements, rotationMatrixAroundAxis(u, Math.PI), tolerance, true)) {
      exactPerpC2Raw.push(u);
    }
    if (checkSymmetryOperation(coords, elements, reflectionMatrixAcrossPlane(u), tolerance)) {
      exactSigmaVRaw.push(u);
    }
  }
  const perpC2s = dedupUnitVectors(exactPerpC2Raw, 0.992);
  const hasSigmaH =
    sigmaList.some((m) => Math.abs(dot3(m, principalAxis!)) > 0.96) ||
    checkSymmetryOperation(coords, elements, reflectionMatrixAcrossPlane(principalAxis), tolerance);
  const sigmaVs = dedupUnitVectors(exactSigmaVRaw, 0.992);
  const hasS2n = checkSymmetryOperation(
    coords,
    elements,
    improperRotationMatrix(principalAxis, Math.PI / principalFold),
    tolerance
  );

  // Build standard frame [frameX, frameY, principalAxis]
  const frameX: Vec3 =
    perpC2s.length > 0
      ? perpC2s[0]
      : sigmaVs.length > 0
        ? sigmaVs[0]
        : normalize3(cross3(principalAxis, Math.abs(principalAxis[0]) < 0.8 ? [1, 0, 0] : [0, 1, 0]));
  const frameY: Vec3 = normalize3(cross3(principalAxis, frameX));

  let pgId = 'C2';
  const hasPerpC2 = perpC2s.length >= principalFold || (principalFold === 2 && perpC2s.length >= 2);

  if (hasPerpC2) {
    // D-family
    if (hasSigmaH) {
      if (principalFold === 6) pgId = 'D6h';
      else if (principalFold === 5) pgId = 'D5h';
      else if (principalFold === 4) pgId = 'D4h';
      else if (principalFold === 3) pgId = 'D3h';
      else pgId = 'D2h';
    } else if (sigmaVs.length >= principalFold) {
      if (principalFold === 5) pgId = 'D5d';
      else if (principalFold === 4) pgId = 'D4d';
      else if (principalFold === 3) pgId = 'D3d';
      else pgId = 'D2d';
    } else {
      if (principalFold === 3) pgId = 'D3';
      else pgId = 'D2';
    }
  } else {
    // C or S family
    if (hasSigmaH) {
      if (principalFold === 3) pgId = 'C3h';
      else if (principalFold === 2) pgId = 'C2h';
      else pgId = 'Cs';
    } else if (sigmaVs.length >= principalFold) {
      if (principalFold === 4) pgId = 'C4v';
      else if (principalFold === 3) pgId = 'C3v';
      else pgId = 'C2v';
    } else if (principalFold === 2 && hasS2n) {
      pgId = 'S4';
    } else {
      pgId = 'C2';
    }
  }

  return buildResult(pgId, centroid, [frameX, frameY, principalAxis], moments, hasInversion, principalFold, tolerance, {
    principalAxis,
    perpC2: perpC2s[0],
    sigmaH: hasSigmaH ? principalAxis : undefined,
    sigmaVs: sigmaVs.slice(0, 3),
  });
}

function buildResult(
  pgId: string,
  centroid: Vec3,
  axes: [Vec3, Vec3, Vec3],
  moments: Vec3,
  hasInversion: boolean,
  principalFold: number,
  tolerance: number,
  detected?: {
    principalAxis?: Vec3;
    perpC2?: Vec3;
    extraAxes?: Vec3[];
    sigmaH?: Vec3;
    sigmaVs?: Vec3[];
  }
): SymmetryDetectionV3 {
  const pg = POINT_GROUPS[pgId] || POINT_GROUPS.C1;
  const [ex, ey, ez] = axes;
  const principalAxes: Mat3 = [
    [ex[0], ey[0], ez[0]],
    [ex[1], ey[1], ez[1]],
    [ex[2], ey[2], ez[2]],
  ];

  // Build aligned 3D symmetry elements in the molecule's real coordinate frame!
  const alignedElements3D: SymmetryElement3D[] = [];
  if (detected?.principalAxis && principalFold >= 2) {
    alignedElements3D.push({
      id: `principal_c${principalFold}`,
      type: 'axis',
      label: principalFold >= 99 ? 'C∞ (Principal Axis)' : `C${principalFold} (Principal Axis)`,
      schoenflies: principalFold >= 99 ? 'C∞' : `C${principalFold}`,
      fold: Math.min(8, principalFold),
      axis: detected.principalAxis,
      color: '#2563EB',
      description: 'Detected principal rotational symmetry axis aligned to the molecular inertia frame',
    });
  }
  if (detected?.perpC2) {
    alignedElements3D.push({
      id: 'perp_c2',
      type: 'axis',
      label: "C2' (Perpendicular Axis)",
      schoenflies: "C2'",
      fold: 2,
      axis: detected.perpC2,
      color: '#0284C7',
      description: 'Perpendicular 2-fold dihedral rotation axis',
    });
  }
  if (detected?.extraAxes) {
    detected.extraAxes.forEach((ax, idx) => {
      alignedElements3D.push({
        id: `extra_axis_${idx}`,
        type: 'axis',
        label: `Symmetry Axis #${idx + 2}`,
        schoenflies: 'Cn',
        fold: 3,
        axis: ax,
        color: '#0284C7',
        description: 'Additional polyhedral/cubic symmetry axis',
      });
    });
  }
  if (detected?.sigmaH) {
    alignedElements3D.push({
      id: 'sigma_h_aligned',
      type: 'plane',
      label: 'σh (Horizontal Mirror)',
      schoenflies: 'σh',
      axis: detected.sigmaH,
      color: '#D97706',
      description: 'Horizontal mirror plane aligned to the molecular coordinates',
    });
  }
  if (detected?.sigmaVs) {
    detected.sigmaVs.forEach((nVec, idx) => {
      alignedElements3D.push({
        id: `sigma_v_${idx}`,
        type: 'plane',
        label: idx === 0 ? 'σv (Vertical Mirror 1)' : `σv' (Vertical Mirror ${idx + 1})`,
        schoenflies: 'σv',
        axis: nVec,
        color: idx % 2 === 0 ? '#0D9488' : '#7C3AED',
        description: 'Vertical / dihedral mirror plane containing the principal axis',
      });
    });
  }
  if (hasInversion) {
    alignedElements3D.push({
      id: 'inversion_center',
      type: 'center',
      label: 'Inversion Center i',
      schoenflies: 'i',
      color: '#DC2626',
      description: 'Exact molecular center of inversion',
    });
  }

  return {
    pointGroupId: pg.id,
    pointGroup: pg,
    centroid,
    principalAxes,
    principalMoments: moments,
    alignedElements3D: alignedElements3D.length > 0 ? alignedElements3D : pg.elements3D,
    hasInversion,
    principalFold,
    rmsdTolerance: tolerance,
  };
}

// ============================================================================
// Exact Wigner Real Spherical-Harmonic Orbital Transformation & Irrep Projection
// ============================================================================

/**
 * Evaluates the 5x5 real d-orbital transformation matrix D^(2)(R) for an orthogonal 3x3 matrix R
 * (including reflections/inversion where det(R) = -1, since d orbitals are gerade).
 */
function buildDOrbitalTransform(R: Mat3): number[][] {
  // Ordered: [dz2, dx2y2, dxy, dxz, dyz]
  const evalD = (x: number, y: number, z: number): number[] => [
    (2 * z * z - x * x - y * y) / Math.sqrt(3),
    x * x - y * y,
    2 * x * y,
    2 * x * z,
    2 * y * z,
  ];
  // Express rotated function f(R^T r) in original basis via exact 26-point octahedral quadrature
  const M: number[][] = Array.from({ length: 5 }, () => new Array<number>(5).fill(0));
  const RT: Mat3 = [
    [R[0][0], R[1][0], R[2][0]],
    [R[0][1], R[1][1], R[2][1]],
    [R[0][2], R[1][2], R[2][2]],
  ];
  const pts: [number, number, number, number][] = [];
  const w1 = 1 / 21;
  const w2 = 4 / 105;
  const w3 = 27 / 840;
  for (const s of [-1, 1]) pts.push([s, 0, 0, w1], [0, s, 0, w1], [0, 0, s, w1]);
  const s2 = Math.SQRT1_2;
  for (const a of [-1, 1])
    for (const b of [-1, 1])
      pts.push([a * s2, b * s2, 0, w2], [a * s2, 0, b * s2, w2], [0, a * s2, b * s2, w2]);
  const s3 = 1 / Math.sqrt(3);
  for (const a of [-1, 1])
    for (const b of [-1, 1])
      for (const c of [-1, 1]) pts.push([a * s3, b * s3, c * s3, w3]);

  // Normalization factor for evalD on S^2: ∫ (evalD_i)^2 dΩ / 4π = 4/15
  const invNorm = 15 / 4;
  for (const [x, y, z, w] of pts) {
    const [rx, ry, rz] = mat3MulVec(RT, [x, y, z]);
    const orig = evalD(x, y, z);
    const rot = evalD(rx, ry, rz);
    for (let i = 0; i < 5; i++) {
      for (let j = 0; j < 5; j++) {
        M[i][j] += w * invNorm * orig[i] * rot[j];
      }
    }
  }
  return M;
}

const D_ORDER: OrbitalBasisType[] = ['dz2', 'dx2y2', 'dxy', 'dxz', 'dyz'];
const P_ORDER: OrbitalBasisType[] = ['px', 'py', 'pz'];

/**
 * Applies a 3D symmetry operation R (in molecular frame around centroid) to an MO coefficient vector `c`
 * and computes its exact character χ(R) = c^T S (D(R) c).
 */
export function computeMOCharacterUnderOperation(
  eht: LowdinEHTResult,
  atoms: AtomNode[],
  centroid: Vec3,
  R: Mat3,
  moIndices: number[],
  tol = 0.28
): number {
  const nAtoms = atoms.length;
  const coords: Vec3[] = atoms.map((a) => [
    a.x3d - centroid[0],
    a.y3d - centroid[1],
    a.z3d - centroid[2],
  ]);

  // Map each atom i -> target atom j = π_R(i) where R * r_i ≈ r_j
  const atomPerm = new Int32Array(nAtoms).fill(-1);
  for (let i = 0; i < nAtoms; i++) {
    const p = mat3MulVec(R, coords[i]);
    let bestJ = -1;
    let bestD = tol;
    for (let j = 0; j < nAtoms; j++) {
      if (atoms[j].element !== atoms[i].element) continue;
      const d = Math.hypot(p[0] - coords[j][0], p[1] - coords[j][1], p[2] - coords[j][2]);
      if (d < bestD) {
        bestD = d;
        bestJ = j;
      }
    }
    if (bestJ < 0) return 0;
    atomPerm[i] = bestJ;
  }

  // Group basis indices by (atomIndex, orbital)
  const basisLookup = new Map<string, number>();
  for (const b of eht.basis) {
    basisLookup.set(`${b.atomIndex}::${b.orbital}`, b.index);
  }

  const Dmat = buildDOrbitalTransform(R);
  const detR =
    R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) -
    R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) +
    R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);

  const nBasis = eht.basis.length;
  let trace = 0;

  for (const k of moIndices) {
    const rc = new Float64Array(nBasis);
    for (let mu = 0; mu < nBasis; mu++) {
      const cMu = eht.coefficients[mu][k];
      if (Math.abs(cMu) < 1e-10) continue;
      const b = eht.basis[mu];
      const targetAtom = atomPerm[b.atomIndex];
      if (b.l === 0) {
        const nu = basisLookup.get(`${targetAtom}::s`);
        if (nu !== undefined) rc[nu] += cMu;
      } else if (b.l === 1) {
        const pIdx = P_ORDER.indexOf(b.orbital);
        for (let q = 0; q < 3; q++) {
          const nu = basisLookup.get(`${targetAtom}::${P_ORDER[q]}`);
          if (nu !== undefined) rc[nu] += R[q][pIdx] * cMu;
        }
      } else if (b.l === 2) {
        const dIdx = D_ORDER.indexOf(b.orbital);
        if (dIdx >= 0) {
          for (let q = 0; q < 5; q++) {
            const nu = basisLookup.get(`${targetAtom}::${D_ORDER[q]}`);
            if (nu !== undefined) rc[nu] += Dmat[q][dIdx] * cMu;
          }
        }
      } else {
        // f-orbital parity approximation under inversion / reflection (ungerade -> * detR)
        const nu = basisLookup.get(`${targetAtom}::${b.orbital}`);
        if (nu !== undefined) rc[nu] += (detR < 0 ? -1 : 1) * cMu;
      }
    }

    // Inner product c_k^T S (R c_k)
    let chi = 0;
    for (let i = 0; i < nBasis; i++) {
      const ci = eht.coefficients[i][k];
      if (Math.abs(ci) < 1e-10) continue;
      let sRc = 0;
      for (let j = 0; j < nBasis; j++) {
        if (rc[j] !== 0) sRc += eht.S[i][j] * rc[j];
      }
      chi += ci * sRc;
    }
    trace += chi;
  }

  return trace;
}

/**
 * Projects every MO (or degenerate MO cluster) from a Löwdin EHT solution onto the exact
 * irreducible representations of the detected 3D point group.
 */
export function projectMOIrrepsExact(
  eht: LowdinEHTResult,
  atoms: AtomNode[],
  sym: SymmetryDetectionV3
): MOIrrepProjection[] {
  const pg = sym.pointGroup;
  const n = eht.energies.length;
  const results: MOIrrepProjection[] = new Array(n);

  // Build representative 3x3 matrix for each symmetry class in the molecule's aligned principal frame
  const ex: Vec3 = [sym.principalAxes[0][0], sym.principalAxes[1][0], sym.principalAxes[2][0]];
  const ey: Vec3 = [sym.principalAxes[0][1], sym.principalAxes[1][1], sym.principalAxes[2][1]];
  const ez: Vec3 = [sym.principalAxes[0][2], sym.principalAxes[1][2], sym.principalAxes[2][2]];

  const classMatrices: Mat3[] = pg.classes.map((cls, idx) => {
    if (cls.operationType === 'E') return IDENTITY_MAT;
    if (cls.operationType === 'i') return INVERSION_MAT;
    if (cls.operationType === 'sigma_h') return reflectionMatrixAcrossPlane(ez);
    if (cls.operationType === 'sigma_v' || cls.operationType === 'sigma_d') {
      const symText = cls.symbol.toLowerCase();
      if (symText.includes('yz') || cls.id.includes('yz') || idx % 2 === 1) {
        return reflectionMatrixAcrossPlane(ex);
      }
      return reflectionMatrixAcrossPlane(ey);
    }
    if (cls.operationType === 'Cn') {
      const m = cls.symbol.match(/C(\d+)/);
      const fold = m ? parseInt(m[1], 10) : sym.principalFold >= 99 ? 4 : sym.principalFold;
      const isPerp = cls.symbol.includes("'") || cls.id.includes('c2p') || cls.id.includes('C2x') || cls.id.includes('C2y');
      if (cls.id.includes('C2x')) return rotationMatrixAroundAxis(ex, Math.PI);
      if (cls.id.includes('C2y')) return rotationMatrixAroundAxis(ey, Math.PI);
      if (isPerp) return rotationMatrixAroundAxis(ex, Math.PI);
      return rotationMatrixAroundAxis(ez, (2 * Math.PI) / Math.max(2, fold));
    }
    if (cls.operationType === 'Sn') {
      const m = cls.symbol.match(/S(\d+)/);
      const fold = m ? parseInt(m[1], 10) : 4;
      return improperRotationMatrix(ez, (2 * Math.PI) / Math.max(2, fold));
    }
    return IDENTITY_MAT;
  });

  // Group nearly degenerate MOs (within 0.045 eV) so 2D E and 3D T subspaces project cleanly
  let idx = 0;
  while (idx < n) {
    let end = idx + 1;
    while (end < n && Math.abs(eht.energies[end] - eht.energies[idx]) < 0.045 && end - idx < 5) {
      end++;
    }
    const cluster = Array.from({ length: end - idx }, (_, j) => idx + j);
    const chars = classMatrices.map((R) =>
      computeMOCharacterUnderOperation(eht, atoms, sym.centroid, R, cluster)
    );

    // Match cluster character vector to the best irrep (or direct sum) via Euclidean distance
    let bestIrrep = pg.irreps[0];
    let bestErr = Infinity;
    for (const ir of pg.irreps) {
      if (cluster.length === 1 && ir.dimension > 1) continue;
      let err = 0;
      const scale = cluster.length / Math.max(1, ir.dimension);
      for (let c = 0; c < pg.classes.length; c++) {
        const diff = chars[c] - ir.characters[c] * scale;
        err += pg.classes[c].count * diff * diff;
      }
      if (err < bestErr) {
        bestErr = err;
        bestIrrep = ir;
      }
    }

    const h = pg.classes.reduce((s, c) => s + c.count, 0) || pg.order || 1;
    const confidence = Math.max(0.2, Math.min(1, 1 - Math.sqrt(bestErr / h) / Math.max(1, cluster.length * 2)));

    for (const moIdx of cluster) {
      results[moIdx] = {
        moIndex: moIdx,
        irrepId: bestIrrep.id,
        irrepSymbol: bestIrrep.symbol,
        confidence,
        computedCharacters: chars.map((v) => Math.round(v * 100) / 100),
      };
    }

    idx = end;
  }

  return results;
}

// ============================================================================
// Optical Selection Rules (UV-Vis Electronic Dipole, Laporte Parity, IR & Raman)
// ============================================================================

function findDipoleIrreps(pg: PointGroupData): {
  x: string[];
  y: string[];
  z: string[];
} {
  const res = { x: [] as string[], y: [] as string[], z: [] as string[] };
  for (const ir of pg.irreps) {
    const joined = ir.linearFunctions.join(' ').toLowerCase();
    // Strip rotations Rx, Ry, Rz before testing translational dipole components x, y, z
    const noRot = joined.replace(/r[xyz]/g, '');
    if (/\bx\b/.test(noRot) || noRot.includes('(x,') || noRot.includes('x, y')) res.x.push(ir.id);
    if (/\by\b/.test(noRot) || noRot.includes(', y') || noRot.includes('x, y')) res.y.push(ir.id);
    if (/\bz\b/.test(noRot) || noRot.includes(', z)')) res.z.push(ir.id);
  }
  return res;
}

/**
 * Evaluates exact group-theoretic electric-dipole selection rules for a transition ψ_i -> ψ_f:
 *   Γ_i ⊗ Γ_μ ⊗ Γ_f ⊃ Γ_tot-sym  <=>  Γ_i ⊗ Γ_f contains Γ_x, Γ_y, or Γ_z.
 */
export function evaluateElectronicTransition(
  pg: PointGroupData,
  fromMo: MolecularOrbital,
  toMo: MolecularOrbital,
  atoms: AtomNode[]
): ElectronicTransitionRule {
  const irFrom =
    pg.irreps.find((ir) => ir.id === fromMo.irrepId || ir.symbol === fromMo.irrepSymbol) ||
    pg.irreps[0];
  const irTo =
    pg.irreps.find((ir) => ir.id === toMo.irrepId || ir.symbol === toMo.irrepSymbol) ||
    pg.irreps[0];

  const productChars = irFrom.characters.map((c, idx) => c * (irTo.characters[idx] ?? 1));
  const decomp = decomposeReducibleRepresentation(pg, productChars);
  const presentIrrepIds = new Set(
    decomp.results.filter((r) => r.coefficient > 0.25).map((r) => r.irrepId)
  );

  const dipoleIrreps = findDipoleIrreps(pg);
  const allowedPolarizations: ('x' | 'y' | 'z')[] = [];
  if (dipoleIrreps.x.some((id) => presentIrrepIds.has(id))) allowedPolarizations.push('x');
  if (dipoleIrreps.y.some((id) => presentIrrepIds.has(id))) allowedPolarizations.push('y');
  if (dipoleIrreps.z.some((id) => presentIrrepIds.has(id))) allowedPolarizations.push('z');

  const isElectricDipoleAllowed = allowedPolarizations.length > 0;

  // Laporte parity rule check for centrosymmetric groups
  const hasInvClass = pg.classes.some((c) => c.operationType === 'i');
  let laporteStatus: 'allowed' | 'forbidden' | 'not-applicable' = 'not-applicable';
  if (hasInvClass) {
    const iIdx = pg.classes.findIndex((c) => c.operationType === 'i');
    const parityFrom = Math.sign(irFrom.characters[iIdx] ?? 0);
    const parityTo = Math.sign(irTo.characters[iIdx] ?? 0);
    if (parityFrom !== 0 && parityTo !== 0) {
      laporteStatus = parityFrom * parityTo < 0 ? 'allowed' : 'forbidden';
    }
  }

  const deltaEv = Math.max(0.01, Math.abs(toMo.energy - fromMo.energy));
  const wavelengthNm = 1239.841984 / deltaEv;
  const wavenumberCm1 = deltaEv * 8065.54429;

  // Approximate transition dipole magnitude |<ψ_i | r | ψ_f>| via Mulliken atomic-center overlap
  const atomCoords = new Map(atoms.map((a) => [a.id, a]));
  let mx = 0,
    my = 0,
    mz = 0;
  for (const t1 of fromMo.lcao) {
    const a = atomCoords.get(t1.atomId);
    if (!a) continue;
    for (const t2 of toMo.lcao) {
      if (t1.atomId === t2.atomId && t1.orbital === t2.orbital) {
        const w = t1.coeff * t2.coeff;
        mx += w * a.x3d;
        my += w * a.y3d;
        mz += w * a.z3d;
      }
    }
  }
  const transitionDipoleStrength = isElectricDipoleAllowed ? Math.hypot(mx, my, mz) : 0;

  let rationale = '';
  if (isElectricDipoleAllowed) {
    rationale = `Direct product ${irFrom.symbol} ⊗ ${irTo.symbol} = ${decomp.formula} spans dipole component(s) ${allowedPolarizations.join(', ')}.`;
    if (laporteStatus === 'allowed') {
      rationale += ' Parity-allowed (g ↔ u, Laporte allowed).';
    }
  } else if (laporteStatus === 'forbidden') {
    rationale = `Laporte forbidden (g → g or u → u): ${irFrom.symbol} ⊗ ${irTo.symbol} = ${decomp.formula} has no ungerade dipole projection.`;
  } else {
    rationale = `Symmetry forbidden: ${irFrom.symbol} ⊗ ${irTo.symbol} = ${decomp.formula} contains no translational dipole irrep (x, y, z).`;
  }

  return {
    fromMoId: fromMo.id,
    fromLabel: fromMo.label,
    fromIrrepId: irFrom.id,
    fromIrrepSymbol: irFrom.symbol,
    toMoId: toMo.id,
    toLabel: toMo.label,
    toIrrepId: irTo.id,
    toIrrepSymbol: irTo.symbol,
    deltaEv,
    wavelengthNm,
    wavenumberCm1,
    directProductFormula: decomp.formula,
    allowedPolarizations,
    isElectricDipoleAllowed,
    laporteStatus,
    transitionDipoleStrength,
    rationale,
  };
}

/**
 * Generates a summary of frontier electronic transitions (top 3 occupied -> bottom 3 virtual)
 * plus IR/Raman vibrational spectroscopic activity for the active point group.
 */
export function computeOpticalSelectionRules(
  pg: PointGroupData,
  orbitals: MolecularOrbital[],
  atoms: AtomNode[]
): {
  frontierTransitions: ElectronicTransitionRule[];
  vibrational: VibrationalActivitySummary;
} {
  const sorted = [...orbitals].sort((a, b) => a.energy - b.energy);
  const occupied = sorted.filter((m) => m.occupancy > 0).slice(-5).reverse();
  const virtual = sorted.filter((m) => m.occupancy < 2).slice(0, 5);

  const frontierTransitions: ElectronicTransitionRule[] = [];
  for (const occ of occupied) {
    for (const virt of virtual) {
      if (occ.id === virt.id || virt.energy - occ.energy <= 0.05) continue;
      frontierTransitions.push(evaluateElectronicTransition(pg, occ, virt, atoms));
    }
  }
  frontierTransitions.sort((a, b) => a.deltaEv - b.deltaEv);

  const isCentrosymmetric = pg.classes.some((c) => c.operationType === 'i');
  const irActiveIrreps = pg.irreps.filter((ir) => ir.irActive).map((ir) => ir.symbol);
  const ramanActiveIrreps = pg.irreps.filter((ir) => ir.ramanActive).map((ir) => ir.symbol);
  const coincidentIrreps = pg.irreps
    .filter((ir) => ir.irActive && ir.ramanActive)
    .map((ir) => ir.symbol);
  const silentIrreps = pg.irreps
    .filter((ir) => !ir.irActive && !ir.ramanActive)
    .map((ir) => ir.symbol);

  return {
    frontierTransitions: frontierTransitions.slice(0, 8),
    vibrational: {
      isCentrosymmetric,
      mutualExclusionApplies: isCentrosymmetric && coincidentIrreps.length === 0,
      irActiveIrreps,
      ramanActiveIrreps,
      coincidentIrreps,
      silentIrreps,
    },
  };
}

export interface VisibleSpectrumColor {
  r: number;
  g: number;
  b: number;
  rgb: string;
  bandName: string;
}

/**
 * Converts a photon wavelength (nm) in the human visible spectrum [380, 750] nm
 * into an approximate perceptual sRGB color (Dan Bruton spectral model with gentle
 * photopic edge roll-off) and spectral band label. Returns null outside 380–750 nm.
 */
export function wavelengthToVisibleColor(wavelengthNm: number): VisibleSpectrumColor | null {
  if (!Number.isFinite(wavelengthNm) || wavelengthNm < 380 || wavelengthNm > 750) {
    return null;
  }

  let r = 0;
  let g = 0;
  let b = 0;

  if (wavelengthNm >= 380 && wavelengthNm < 440) {
    r = -(wavelengthNm - 440) / (440 - 380);
    g = 0;
    b = 1;
  } else if (wavelengthNm >= 440 && wavelengthNm < 490) {
    r = 0;
    g = (wavelengthNm - 440) / (490 - 440);
    b = 1;
  } else if (wavelengthNm >= 490 && wavelengthNm < 510) {
    r = 0;
    g = 1;
    b = -(wavelengthNm - 510) / (510 - 490);
  } else if (wavelengthNm >= 510 && wavelengthNm < 580) {
    r = (wavelengthNm - 510) / (580 - 510);
    g = 1;
    b = 0;
  } else if (wavelengthNm >= 580 && wavelengthNm < 645) {
    r = 1;
    g = -(wavelengthNm - 645) / (645 - 580);
    b = 0;
  } else {
    r = 1;
    g = 0;
    b = 0;
  }

  // Gentle photopic roll-off at deep violet & deep red ends while keeping hues distinct
  let intensity = 1.0;
  if (wavelengthNm < 420) {
    intensity = 0.48 + 0.52 * ((wavelengthNm - 380) / (420 - 380));
  } else if (wavelengthNm > 680) {
    intensity = 0.48 + 0.52 * ((750 - wavelengthNm) / (750 - 680));
  }

  const gamma = 0.8;
  const R = r > 0 ? Math.round(255 * Math.pow(r * intensity, gamma)) : 0;
  const G = g > 0 ? Math.round(255 * Math.pow(g * intensity, gamma)) : 0;
  const B = b > 0 ? Math.round(255 * Math.pow(b * intensity, gamma)) : 0;

  let bandName = 'Red';
  if (wavelengthNm < 435) bandName = 'Violet';
  else if (wavelengthNm < 485) bandName = 'Blue';
  else if (wavelengthNm < 510) bandName = 'Cyan';
  else if (wavelengthNm < 565) bandName = 'Green';
  else if (wavelengthNm < 590) bandName = 'Yellow';
  else if (wavelengthNm < 625) bandName = 'Orange';

  return {
    r: R,
    g: G,
    b: B,
    rgb: `rgb(${R}, ${G}, ${B})`,
    bandName,
  };
}

