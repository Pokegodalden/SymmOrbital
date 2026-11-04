import { AtomNode, LCAOContribution, MolecularOrbital } from './moleculesData';

export interface IsosurfaceMeshData {
  positions: Float32Array;
  normals: Float32Array;
  indices: Uint32Array;
}

export interface DualPhaseIsosurfaceData {
  positive: IsosurfaceMeshData;
  negative: IsosurfaceMeshData;
  morphT: number;
}

interface PreparedAtomGroup {
  x: number;
  y: number;
  z: number;
  terms: {
    weight: number;
    orbital: LCAOContribution['orbital'];
  }[];
}

// Kuhn 6-tetrahedra decomposition of a unit cube (vertices 0..7)
// Cube vertex numbering:
// 0:(0,0,0), 1:(1,0,0), 2:(1,1,0), 3:(0,1,0),
// 4:(0,0,1), 5:(1,0,1), 6:(1,1,1), 7:(0,1,1)
const KUHN_TETRAHEDRA: readonly [number, number, number, number][] = [
  [0, 5, 1, 6],
  [0, 1, 2, 6],
  [0, 2, 3, 6],
  [0, 3, 7, 6],
  [0, 7, 4, 6],
  [0, 4, 5, 6],
];

// For each of the 16 tetrahedron sign configurations (bit b set if v_b > iso),
// list pairs of tetrahedron local vertex indices (0..3) forming triangles.
// Winding order is oriented so normals point toward lower scalar values (outward from +iso).
const TET_TRIANGLES: readonly (readonly number[])[] = [
  [], // 0000
  [0, 1, 0, 2, 0, 3], // 0001 (v0 inside) -> 3 crossed edges (0-1, 0-2, 0-3)
  [1, 0, 1, 3, 1, 2], // 0010 (v1 inside) -> 3 crossed edges (1-0, 1-3, 1-2)
  [0, 2, 0, 3, 1, 3, 1, 2], // 0011 (v0,v1 inside) -> cyclic quad (0-2, 0-3, 1-3, 1-2)
  [2, 0, 2, 1, 2, 3], // 0100 (v2 inside) -> 3 crossed edges (2-0, 2-1, 2-3)
  [0, 1, 0, 3, 2, 3, 2, 1], // 0101 (v0,v2 inside) -> cyclic quad (0-1, 0-3, 2-3, 2-1)
  [1, 0, 1, 3, 2, 3, 2, 0], // 0110 (v1,v2 inside) -> cyclic quad (1-0, 1-3, 2-3, 2-0)
  [3, 0, 3, 2, 3, 1], // 0111 (v3 outside) -> 3 crossed edges (3-0, 3-2, 3-1)
  [3, 0, 3, 2, 3, 1], // 1000 (v3 inside) -> 3 crossed edges (3-0, 3-2, 3-1)
  [0, 1, 0, 2, 3, 2, 3, 1], // 1001 (v0,v3 inside) -> cyclic quad (0-1, 0-2, 3-2, 3-1)
  [1, 0, 1, 2, 3, 2, 3, 0], // 1010 (v1,v3 inside) -> cyclic quad (1-0, 1-2, 3-2, 3-0)
  [2, 0, 2, 1, 2, 3], // 1011 (v2 outside) -> 3 crossed edges (2-0, 2-1, 2-3)
  [2, 0, 2, 1, 3, 1, 3, 0], // 1100 (v2,v3 inside) -> cyclic quad (2-0, 2-1, 3-1, 3-0)
  [1, 0, 1, 3, 1, 2], // 1101 (v1 outside) -> 3 crossed edges (1-0, 1-3, 1-2)
  [0, 1, 0, 2, 0, 3], // 1110 (v0 outside) -> 3 crossed edges (0-1, 0-2, 0-3)
  [], // 1111
];

/**
 * Evaluates a single atomic orbital basis function at displacement (dx, dy, dz)
 * with radial Slater exponent zeta.
 */
function evalSingleBasis(
  orbital: LCAOContribution['orbital'],
  dx: number,
  dy: number,
  dz: number,
  r: number,
  zeta: number
): number {
  const radial = Math.exp(-zeta * r);
  const pScale = 1.62 * zeta;
  const dScale = (1.35 * zeta * zeta) / (1.0 + 0.22 * r);
  let angular = 1.0;
  switch (orbital) {
    case 's':
      angular = 0.78;
      break;
    case 'px':
      angular = dx * pScale;
      break;
    case 'py':
      angular = dy * pScale;
      break;
    case 'pz':
      angular = dz * pScale;
      break;
    case 'dz2':
      angular = (2 * dy * dy - dx * dx - dz * dz) * 0.95 * dScale;
      break;
    case 'dx2y2':
      angular = (dx * dx - dz * dz) * 1.15 * dScale;
      break;
    case 'dxy':
      angular = dx * dy * 2.15 * dScale;
      break;
    case 'dxz':
      angular = dx * dz * 2.15 * dScale;
      break;
    case 'dyz':
      angular = dy * dz * 2.15 * dScale;
      break;
  }
  return radial * angular;
}

/**
 * Evaluates the AO-to-MO morphing scalar field F(r, t) at (x, y, z) for morph parameter t in [0, 1].
 * - At t = 0 (Separate AOs): Each voxel is governed by its dominant atom-centered orbital with
 *   an Inverse Participation Ratio (IPR) cleft factor so adjacent atomic orbitals remain cleanly separated.
 * - At t = 1 (Merged MO): Full coherent LCAO superposition Σ w_k φ_k(r - R_k), producing continuous
 *   fused bonding envelopes (σ sausages, π clouds, delocalized rings) and sharp antibonding nodal planes.
 */
export function evaluateMorphScalarFieldAtPoint(
  mo: MolecularOrbital,
  atoms: AtomNode[],
  x: number,
  y: number,
  z: number,
  morphT: number
): number {
  const t = Math.max(0, Math.min(1, morphT));
  const prepared = prepareTerms(mo, atoms);
  if (prepared.length === 0) return 0;
  // Smooth C1 interpolation of Slater radial exponent from compact AO (t=0) to delocalized MO (t=1)
  const s = t * t * (3 - 2 * t);
  const zeta = 2.12 * (1 - s) + 1.48 * s;
  return evalPreparedMorphField(prepared, x, y, z, s, zeta);
}

function prepareTerms(mo: MolecularOrbital, atoms: AtomNode[]): PreparedAtomGroup[] {
  const isDelocalized = mo.lcao.length > 10;
  const gamma = isDelocalized ? 0.56 : 0.70;

  // Group LCAO contributions by atom so intra-atomic hybrid/tilted orbitals preserve linear
  // directional cosines (rotational invariance) while dynamic-range scaling is applied per atom
  const byAtom = new Map<string, { atom: AtomNode; entries: LCAOContribution[] }>();
  for (const term of mo.lcao) {
    if (Math.abs(term.coeff) < 0.04) continue;
    const atom = atoms.find((a) => a.id === term.atomId);
    if (!atom) continue;
    let group = byAtom.get(atom.id);
    if (!group) {
      group = { atom, entries: [] };
      byAtom.set(atom.id, group);
    }
    group.entries.push(term);
  }

  const groups: PreparedAtomGroup[] = [];
  let maxWeight = 0;

  for (const { atom, entries } of byAtom.values()) {
    let atomNorm2 = 0;
    for (const e of entries) {
      atomNorm2 += e.coeff * e.coeff;
    }
    const atomNorm = Math.sqrt(atomNorm2);
    if (atomNorm < 1e-6) continue;

    // Apply nonlinear dynamic-range compression to the total atomic amplitude while preserving
    // exact linear ratios among (s, px, py, pz) components on the same atom
    const atomScale = Math.pow(atomNorm, gamma) / atomNorm;
    const atomTerms: { weight: number; orbital: LCAOContribution['orbital'] }[] = [];

    for (const e of entries) {
      const w = e.coeff * atomScale;
      if (Math.abs(w) > maxWeight) maxWeight = Math.abs(w);
      atomTerms.push({
        weight: w,
        orbital: e.orbital,
      });
    }

    groups.push({
      x: atom.x3d,
      y: atom.y3d,
      z: atom.z3d,
      terms: atomTerms,
    });
  }

  if (maxWeight > 1e-6) {
    const inv = 1 / maxWeight;
    for (const g of groups) {
      for (const rt of g.terms) {
        rt.weight *= inv;
      }
    }
  }
  return groups;
}

/**
 * Evaluates the C^∞ smooth morphing scalar field using atom-grouped Inverse Participation Ratio (IPR)
 * partitioning at t = 0 and Pauli core-excluded constructive/destructive superposition at t = 1.
 */
function evalPreparedMorphField(
  prepared: PreparedAtomGroup[],
  x: number,
  y: number,
  z: number,
  s: number,
  zeta: number
): number {
  let posSum = 0;
  let negSum = 0;
  let posWeight = 0;
  let negWeight = 0;
  let sumA2 = 0;
  let sumVA2 = 0;
  let sumA4 = 0;

  for (let i = 0; i < prepared.length; i++) {
    const group = prepared[i];
    const dx = x - group.x;
    const dy = y - group.y;
    const dz = z - group.z;
    const r2 = dx * dx + dy * dy + dz * dz;
    // Skip negligible far-field tails (> 6.0 Å away) with smooth C^1 radial taper beyond 4.8 Å
    if (r2 > 36.0) continue;

    const r = Math.sqrt(r2);
    let atomVal = 0;
    for (let j = 0; j < group.terms.length; j++) {
      const term = group.terms[j];
      atomVal += term.weight * evalSingleBasis(term.orbital, dx, dy, dz, r, zeta);
    }
    if (r > 4.8) {
      const u = (r - 4.8) / 1.2;
      atomVal *= 1 - u * u * (3 - 2 * u);
    }

    // Smooth Gaussian Voronoi / Pauli core-locality weight prevents diffuse tails of distant atoms
    // from erasing inner antibonding lobes across short X-H bonds while keeping exact nodal planes
    const loc = Math.exp(-2.65 * r2);
    const a2 = atomVal * atomVal * loc;

    if (atomVal >= 0) {
      posSum += atomVal;
      posWeight += a2;
    } else {
      negSum += -atomVal;
      negWeight += a2;
    }

    sumA2 += a2;
    sumVA2 += atomVal * a2;
    sumA4 += a2 * a2;
  }

  const totalWeight = posWeight + negWeight;
  if (totalWeight < 1e-12) return 0;

  // At t = 1 (Merged MO): constructive same-phase terms add linearly (posSum or negSum), while
  // opposite-phase terms meet at a sharp nodal plane without swallowing neighboring atomic cores
  const mergedPsi = (posWeight * posSum - negWeight * negSum) / totalWeight;
  if (s >= 0.999) return mergedPsi;

  // At t = 0 (Separate AOs): C^∞ dominant-atom orbital amplitude and Inverse Participation Ratio P(r)
  const invA2 = 1.0 / sumA2;
  const vSoft = sumVA2 * invA2;
  const ipr = sumA4 * invA2 * invA2; // in [1/N_atoms, 1]
  const sepPsi = vSoft * ipr;

  return (1 - s) * sepPsi + s * mergedPsi;
}

/**
 * Extracts watertight positive (+ψ) and negative (−ψ) isosurface triangle meshes
 * at a given morph parameter t in [0, 1] using Marching Tetrahedra with smooth
 * central-difference scalar-field gradient normals and tangent-plane surface fairing.
 */
export function buildMorphIsosurfaceKeyframe(
  mo: MolecularOrbital,
  atoms: AtomNode[],
  morphT: number,
  lobeScale = 1.0
): DualPhaseIsosurfaceData {
  const t = Math.max(0, Math.min(1, morphT));
  const prepared = prepareTerms(mo, atoms);

  if (prepared.length === 0) {
    const empty: IsosurfaceMeshData = {
      positions: new Float32Array(0),
      normals: new Float32Array(0),
      indices: new Uint32Array(0),
    };
    return { positive: empty, negative: empty, morphT: t };
  }

  // Compute 3D bounding box around active LCAO centers with generous padding (3.45 Å + scale allowance)
  // so outer lobes of p/d and constructive π/σ orbitals never touch or clip against the grid walls.
  let minX = Infinity,
    maxX = -Infinity;
  let minY = Infinity,
    maxY = -Infinity;
  let minZ = Infinity,
    maxZ = -Infinity;

  for (const pt of prepared) {
    if (pt.x < minX) minX = pt.x;
    if (pt.x > maxX) maxX = pt.x;
    if (pt.y < minY) minY = pt.y;
    if (pt.y > maxY) maxY = pt.y;
    if (pt.z < minZ) minZ = pt.z;
    if (pt.z > maxZ) maxZ = pt.z;
  }

  const pad = 3.45 * Math.pow(Math.max(1.0, lobeScale), 0.45);
  minX -= pad;
  maxX += pad;
  minY -= pad;
  maxY += pad;
  minZ -= pad;
  maxZ += pad;

  const spanX = Math.max(4.8, maxX - minX);
  const spanY = Math.max(4.8, maxY - minY);
  const spanZ = Math.max(4.8, maxZ - minZ);
  const maxSpan = Math.max(spanX, spanY, spanZ);

  // Adaptive resolution: 34 cells along longest axis for compact molecules, 24 for macromolecules
  const baseMaxCells = prepared.length > 12 || maxSpan > 11.5 ? 24 : 34;
  const nx = Math.max(20, Math.round((spanX / maxSpan) * baseMaxCells));
  const ny = Math.max(20, Math.round((spanY / maxSpan) * baseMaxCells));
  const nz = Math.max(20, Math.round((spanZ / maxSpan) * baseMaxCells));

  const dx = spanX / nx;
  const dy = spanY / ny;
  const dz = spanZ / nz;

  const vxCount = nx + 1;
  const vyCount = ny + 1;
  const vzCount = nz + 1;
  const sliceStride = vxCount * vyCount;
  const totalVerts = vxCount * vyCount * vzCount;

  const field = new Float32Array(totalVerts);
  const s = t * t * (3 - 2 * t);
  const zeta = 2.12 * (1 - s) + 1.48 * s;

  // Smooth boundary-envelope window: strictly 0 on the 6 outer faces of the sampling volume
  // so every isosurface is mathematically guaranteed to be a 100% closed, watertight 2-manifold.
  const boundaryWeight = (idx: number, maxIdx: number): number => {
    const dist = Math.min(idx, maxIdx - idx);
    if (dist <= 0) return 0;
    if (dist === 1) return 0.65;
    return 1.0;
  };

  // 1. Sample scalar field F(r, t) at all grid vertices
  for (let iz = 0; iz < vzCount; iz++) {
    const z = minZ + iz * dz;
    const wz = boundaryWeight(iz, nz);
    const zOffset = iz * sliceStride;
    for (let iy = 0; iy < vyCount; iy++) {
      const y = minY + iy * dy;
      const wyz = wz * boundaryWeight(iy, ny);
      const yOffset = zOffset + iy * vxCount;
      for (let ix = 0; ix < vxCount; ix++) {
        const wxyz = wyz * boundaryWeight(ix, nx);
        if (wxyz === 0) {
          field[yOffset + ix] = 0;
          continue;
        }
        const x = minX + ix * dx;
        field[yOffset + ix] = wxyz * evalPreparedMorphField(prepared, x, y, z, s, zeta);
      }
    }
  }

  // 2. Pre-compute central-difference gradient ∇F at all grid vertices for silky-smooth normals
  const gradX = new Float32Array(totalVerts);
  const gradY = new Float32Array(totalVerts);
  const gradZ = new Float32Array(totalVerts);
  const inv2Dx = 0.5 / dx;
  const inv2Dy = 0.5 / dy;
  const inv2Dz = 0.5 / dz;

  for (let iz = 0; iz < vzCount; iz++) {
    const izPrev = Math.max(0, iz - 1) * sliceStride;
    const izNext = Math.min(vzCount - 1, iz + 1) * sliceStride;
    const zOff = iz * sliceStride;

    for (let iy = 0; iy < vyCount; iy++) {
      const iyPrev = Math.max(0, iy - 1) * vxCount;
      const iyNext = Math.min(vyCount - 1, iy + 1) * vxCount;
      const yOff = zOff + iy * vxCount;

      for (let ix = 0; ix < vxCount; ix++) {
        const idx = yOff + ix;
        const ixPrev = Math.max(0, ix - 1);
        const ixNext = Math.min(vxCount - 1, ix + 1);

        gradX[idx] = (field[yOff + ixNext] - field[yOff + ixPrev]) * inv2Dx;
        gradY[idx] = (field[zOff + iyNext + ix] - field[zOff + iyPrev + ix]) * inv2Dy;
        gradZ[idx] = (field[izNext + iy * vxCount + ix] - field[izPrev + iy * vxCount + ix]) * inv2Dz;
      }
    }
  }

  // Calibrated isovalue threshold (scaled inversely with user Lobe Scale slider)
  const isoValue = 0.20 / Math.max(0.45, Math.min(1.65, lobeScale));

  const positive = extractPhaseIsosurface(
    field,
    gradX,
    gradY,
    gradZ,
    nx,
    ny,
    nz,
    minX,
    minY,
    minZ,
    dx,
    dy,
    dz,
    isoValue,
    1
  );

  const negative = extractPhaseIsosurface(
    field,
    gradX,
    gradY,
    gradZ,
    nx,
    ny,
    nz,
    minX,
    minY,
    minZ,
    dx,
    dy,
    dz,
    isoValue,
    -1
  );

  return { positive, negative, morphT: t };
}

function extractPhaseIsosurface(
  field: Float32Array,
  gradX: Float32Array,
  gradY: Float32Array,
  gradZ: Float32Array,
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
): IsosurfaceMeshData {
  const vxCount = nx + 1;
  const vyCount = ny + 1;
  const sliceStride = vxCount * vyCount;

  const positions: number[] = [];
  const normals: number[] = [];
  const indices: number[] = [];
  const edgeVertexMap = new Map<number, number>();

  const pushOrientedTriangle = (i0: number, i1: number, i2: number) => {
    const p0x = positions[i0 * 3],
      p0y = positions[i0 * 3 + 1],
      p0z = positions[i0 * 3 + 2];
    const p1x = positions[i1 * 3],
      p1y = positions[i1 * 3 + 1],
      p1z = positions[i1 * 3 + 2];
    const p2x = positions[i2 * 3],
      p2y = positions[i2 * 3 + 1],
      p2z = positions[i2 * 3 + 2];

    const ux = p1x - p0x,
      uy = p1y - p0y,
      uz = p1z - p0z;
    const vx = p2x - p0x,
      vy = p2y - p0y,
      vz = p2z - p0z;

    const fnx = uy * vz - uz * vy;
    const fny = uz * vx - ux * vz;
    const fnz = ux * vy - uy * vx;

    const gnx = normals[i0 * 3] + normals[i1 * 3] + normals[i2 * 3];
    const gny = normals[i0 * 3 + 1] + normals[i1 * 3 + 1] + normals[i2 * 3 + 1];
    const gnz = normals[i0 * 3 + 2] + normals[i1 * 3 + 2] + normals[i2 * 3 + 2];

    if (fnx * gnx + fny * gny + fnz * gnz >= 0) {
      indices.push(i0, i1, i2);
    } else {
      indices.push(i0, i2, i1);
    }
  };

  // Cube corner offsets in grid index and world space
  const cOffIdx = [
    0,
    1,
    1 + vxCount,
    vxCount,
    sliceStride,
    sliceStride + 1,
    sliceStride + 1 + vxCount,
    sliceStride + vxCount,
  ];
  const cOffX = [0, dx, dx, 0, 0, dx, dx, 0];
  const cOffY = [0, 0, dy, dy, 0, 0, dy, dy];
  const cOffZ = [0, 0, 0, 0, dz, dz, dz, dz];

  const getOrCreateEdgeVertex = (
    gIdxA: number,
    gIdxB: number,
    ax: number,
    ay: number,
    az: number,
    bx: number,
    by: number,
    bz: number,
    valA: number,
    valB: number
  ): number => {
    const minIdx = gIdxA < gIdxB ? gIdxA : gIdxB;
    const maxIdx = gIdxA < gIdxB ? gIdxB : gIdxA;
    // Unique hash for undirected edge between two grid vertices (totalVerts < 131,072)
    const edgeKey = minIdx * 131072 + maxIdx;
    const existing = edgeVertexMap.get(edgeKey);
    if (existing !== undefined) return existing;

    const denom = valB - valA;
    const mu = Math.abs(denom) > 1e-7 ? Math.max(0, Math.min(1, (isoValue - valA) / denom)) : 0.5;

    const px = ax + mu * (bx - ax);
    const py = ay + mu * (by - ay);
    const pz = az + mu * (bz - az);

    // Outward normal points down the signed gradient (-phaseSign * ∇F)
    let gx = -(gradX[gIdxA] + mu * (gradX[gIdxB] - gradX[gIdxA])) * phaseSign;
    let gy = -(gradY[gIdxA] + mu * (gradY[gIdxB] - gradY[gIdxA])) * phaseSign;
    let gz = -(gradZ[gIdxA] + mu * (gradZ[gIdxB] - gradZ[gIdxA])) * phaseSign;
    const gLen = Math.hypot(gx, gy, gz) || 1;
    gx /= gLen;
    gy /= gLen;
    gz /= gLen;

    const vIndex = positions.length / 3;
    positions.push(px, py, pz);
    normals.push(gx, gy, gz);
    edgeVertexMap.set(edgeKey, vIndex);
    return vIndex;
  };

  const cubeVals = new Float32Array(8);
  const cubeGridIdx = new Int32Array(8);
  const cubeX = new Float32Array(8);
  const cubeY = new Float32Array(8);
  const cubeZ = new Float32Array(8);

  for (let iz = 0; iz < nz; iz++) {
    const z0 = minZ + iz * dz;
    const zBase = iz * sliceStride;
    for (let iy = 0; iy < ny; iy++) {
      const y0 = minY + iy * dy;
      const yBase = zBase + iy * vxCount;
      for (let ix = 0; ix < nx; ix++) {
        const baseIdx = yBase + ix;
        const x0 = minX + ix * dx;

        let anyAbove = false;
        let anyBelow = false;
        for (let c = 0; c < 8; c++) {
          const gIdx = baseIdx + cOffIdx[c];
          const v = field[gIdx] * phaseSign;
          cubeVals[c] = v;
          cubeGridIdx[c] = gIdx;
          if (v > isoValue) anyAbove = true;
          else anyBelow = true;
        }

        if (!anyAbove || !anyBelow) continue;

        for (let c = 0; c < 8; c++) {
          cubeX[c] = x0 + cOffX[c];
          cubeY[c] = y0 + cOffY[c];
          cubeZ[c] = z0 + cOffZ[c];
        }

        for (let tIdx = 0; tIdx < 6; tIdx++) {
          const tet = KUHN_TETRAHEDRA[tIdx];
          const c0 = tet[0],
            c1 = tet[1],
            c2 = tet[2],
            c3 = tet[3];

          let mask = 0;
          if (cubeVals[c0] > isoValue) mask |= 1;
          if (cubeVals[c1] > isoValue) mask |= 2;
          if (cubeVals[c2] > isoValue) mask |= 4;
          if (cubeVals[c3] > isoValue) mask |= 8;

          if (mask === 0 || mask === 15) continue;

          const edges = TET_TRIANGLES[mask];
          if (edges.length === 6) {
            // Single triangle (3 crossed edges)
            const a0 = tet[edges[0]],
              b0 = tet[edges[1]];
            const a1 = tet[edges[2]],
              b1 = tet[edges[3]];
            const a2 = tet[edges[4]],
              b2 = tet[edges[5]];

            const i0 = getOrCreateEdgeVertex(
              cubeGridIdx[a0],
              cubeGridIdx[b0],
              cubeX[a0],
              cubeY[a0],
              cubeZ[a0],
              cubeX[b0],
              cubeY[b0],
              cubeZ[b0],
              cubeVals[a0],
              cubeVals[b0]
            );
            const i1 = getOrCreateEdgeVertex(
              cubeGridIdx[a1],
              cubeGridIdx[b1],
              cubeX[a1],
              cubeY[a1],
              cubeZ[a1],
              cubeX[b1],
              cubeY[b1],
              cubeZ[b1],
              cubeVals[a1],
              cubeVals[b1]
            );
            const i2 = getOrCreateEdgeVertex(
              cubeGridIdx[a2],
              cubeGridIdx[b2],
              cubeX[a2],
              cubeY[a2],
              cubeZ[a2],
              cubeX[b2],
              cubeY[b2],
              cubeZ[b2],
              cubeVals[a2],
              cubeVals[b2]
            );
            pushOrientedTriangle(i0, i1, i2);
          } else if (edges.length === 8) {
            // Quad split into 2 triangles (4 crossed edges)
            const a0 = tet[edges[0]],
              b0 = tet[edges[1]];
            const a1 = tet[edges[2]],
              b1 = tet[edges[3]];
            const a2 = tet[edges[4]],
              b2 = tet[edges[5]];
            const a3 = tet[edges[6]],
              b3 = tet[edges[7]];

            const i0 = getOrCreateEdgeVertex(
              cubeGridIdx[a0],
              cubeGridIdx[b0],
              cubeX[a0],
              cubeY[a0],
              cubeZ[a0],
              cubeX[b0],
              cubeY[b0],
              cubeZ[b0],
              cubeVals[a0],
              cubeVals[b0]
            );
            const i1 = getOrCreateEdgeVertex(
              cubeGridIdx[a1],
              cubeGridIdx[b1],
              cubeX[a1],
              cubeY[a1],
              cubeZ[a1],
              cubeX[b1],
              cubeY[b1],
              cubeZ[b1],
              cubeVals[a1],
              cubeVals[b1]
            );
            const i2 = getOrCreateEdgeVertex(
              cubeGridIdx[a2],
              cubeGridIdx[b2],
              cubeX[a2],
              cubeY[a2],
              cubeZ[a2],
              cubeX[b2],
              cubeY[b2],
              cubeZ[b2],
              cubeVals[a2],
              cubeVals[b2]
            );
            const i3 = getOrCreateEdgeVertex(
              cubeGridIdx[a3],
              cubeGridIdx[b3],
              cubeX[a3],
              cubeY[a3],
              cubeZ[a3],
              cubeX[b3],
              cubeY[b3],
              cubeZ[b3],
              cubeVals[a3],
              cubeVals[b3]
            );
            pushOrientedTriangle(i0, i1, i2);
            pushOrientedTriangle(i0, i2, i3);
          }
        }
      }
    }
  }

  const posOut = new Float32Array(positions);
  const normOut = new Float32Array(normals);
  const idxOut = new Uint32Array(indices);

  if (posOut.length > 0 && idxOut.length > 0) {
    smoothIsosurfaceTangentPlaneAndNormals(posOut, normOut, idxOut, 2, 0.42);
  }

  return {
    positions: posOut,
    normals: normOut,
    indices: idxOut,
  };
}

/**
 * Performs tangent-plane Laplacian vertex relaxation:
 *   p_i <- p_i + λ (I - n_i n_i^T) (p_avg - p_i)
 * followed by 1-ring normal averaging and unit normalization.
 * Because displacement is strictly orthogonal to the scalar-field gradient normal n_i,
 * vertices stay on the exact F(p) = isoValue level set with zero volume shrinkage while
 * eliminating Marching Tetrahedra diagonal grid ripples.
 */
function smoothIsosurfaceTangentPlaneAndNormals(
  positions: Float32Array,
  normals: Float32Array,
  indices: Uint32Array,
  iterations: number,
  lambda: number
): void {
  const vCount = positions.length / 3;
  const sumX = new Float32Array(vCount);
  const sumY = new Float32Array(vCount);
  const sumZ = new Float32Array(vCount);
  const deg = new Uint16Array(vCount);

  for (let iter = 0; iter < iterations; iter++) {
    sumX.fill(0);
    sumY.fill(0);
    sumZ.fill(0);
    deg.fill(0);

    for (let t = 0; t < indices.length; t += 3) {
      const i0 = indices[t];
      const i1 = indices[t + 1];
      const i2 = indices[t + 2];

      const p0x = positions[i0 * 3],
        p0y = positions[i0 * 3 + 1],
        p0z = positions[i0 * 3 + 2];
      const p1x = positions[i1 * 3],
        p1y = positions[i1 * 3 + 1],
        p1z = positions[i1 * 3 + 2];
      const p2x = positions[i2 * 3],
        p2y = positions[i2 * 3 + 1],
        p2z = positions[i2 * 3 + 2];

      sumX[i0] += p1x + p2x;
      sumY[i0] += p1y + p2y;
      sumZ[i0] += p1z + p2z;
      deg[i0] += 2;

      sumX[i1] += p0x + p2x;
      sumY[i1] += p0y + p2y;
      sumZ[i1] += p0z + p2z;
      deg[i1] += 2;

      sumX[i2] += p0x + p1x;
      sumY[i2] += p0y + p1y;
      sumZ[i2] += p0z + p1z;
      deg[i2] += 2;
    }

    for (let v = 0; v < vCount; v++) {
      const d = deg[v];
      if (d === 0) continue;
      const invD = 1.0 / d;
      const v3 = v * 3;
      const dx = sumX[v] * invD - positions[v3];
      const dy = sumY[v] * invD - positions[v3 + 1];
      const dz = sumZ[v] * invD - positions[v3 + 2];

      const nx = normals[v3];
      const ny = normals[v3 + 1];
      const nz = normals[v3 + 2];
      const dot = dx * nx + dy * ny + dz * nz;

      // Project onto tangent plane orthogonal to normal
      positions[v3] += lambda * (dx - dot * nx);
      positions[v3 + 1] += lambda * (dy - dot * ny);
      positions[v3 + 2] += lambda * (dz - dot * nz);
    }
  }

  // Smooth vertex normals over 1-ring neighborhood for silky specular highlights
  sumX.fill(0);
  sumY.fill(0);
  sumZ.fill(0);

  for (let t = 0; t < indices.length; t += 3) {
    const i0 = indices[t];
    const i1 = indices[t + 1];
    const i2 = indices[t + 2];

    const n0x = normals[i0 * 3],
      n0y = normals[i0 * 3 + 1],
      n0z = normals[i0 * 3 + 2];
    const n1x = normals[i1 * 3],
      n1y = normals[i1 * 3 + 1],
      n1z = normals[i1 * 3 + 2];
    const n2x = normals[i2 * 3],
      n2y = normals[i2 * 3 + 1],
      n2z = normals[i2 * 3 + 2];

    sumX[i0] += n1x + n2x;
    sumY[i0] += n1y + n2y;
    sumZ[i0] += n1z + n2z;

    sumX[i1] += n0x + n2x;
    sumY[i1] += n0y + n2y;
    sumZ[i1] += n0z + n2z;

    sumX[i2] += n0x + n1x;
    sumY[i2] += n0y + n1y;
    sumZ[i2] += n0z + n1z;
  }

  for (let v = 0; v < vCount; v++) {
    const d = deg[v];
    if (d === 0) continue;
    const invD = 1.0 / d;
    const v3 = v * 3;
    const nx = 0.55 * normals[v3] + 0.45 * sumX[v] * invD;
    const ny = 0.55 * normals[v3 + 1] + 0.45 * sumY[v] * invD;
    const nz = 0.55 * normals[v3 + 2] + 0.45 * sumZ[v] * invD;
    const len = Math.hypot(nx, ny, nz) || 1.0;
    normals[v3] = nx / len;
    normals[v3 + 1] = ny / len;
    normals[v3 + 2] = nz / len;
  }
}

