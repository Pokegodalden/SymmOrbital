import {
  AtomNode,
  BondEdge,
  ELEMENT_PROPERTIES,
  ElementSymbol,
} from './moleculesData';
import { ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM } from './quantumEngineV3';

export interface ParsedXYZResult {
  title: string;
  atoms: AtomNode[];
  bonds: BondEdge[];
}

export interface GaussianCubeData {
  title: string;
  comment: string;
  atoms: AtomNode[];
  bonds: BondEdge[];
  originAngstrom: [number, number, number];
  nx: number;
  ny: number;
  nz: number;
  vxAngstrom: [number, number, number];
  vyAngstrom: [number, number, number];
  vzAngstrom: [number, number, number];
  values: Float32Array;
  minVal: number;
  maxVal: number;
}

const ATOMIC_NUMBER_TO_SYMBOL: Record<number, ElementSymbol> = (() => {
  const map: Record<number, ElementSymbol> = {};
  for (const [sym, props] of Object.entries(ELEMENT_PROPERTIES)) {
    map[props.atomicNumber] = sym as ElementSymbol;
  }
  return map;
})();

function normalizeElementSymbol(raw: string): ElementSymbol | null {
  const trimmed = raw.trim();
  if (!trimmed) return null;
  if (/^\d+$/.test(trimmed)) {
    const z = parseInt(trimmed, 10);
    return ATOMIC_NUMBER_TO_SYMBOL[z] || null;
  }
  const normalized =
    trimmed.charAt(0).toUpperCase() + trimmed.slice(1).toLowerCase();
  if (normalized in ELEMENT_PROPERTIES) {
    return normalized as ElementSymbol;
  }
  return null;
}

/**
 * Infers covalent and coordination bonds from 3D Cartesian coordinates using covalent radii.
 */
export function perceiveBondsFrom3D(atoms: AtomNode[]): BondEdge[] {
  const bonds: BondEdge[] = [];
  let idx = 1;
  for (let i = 0; i < atoms.length; i++) {
    const p1 = ELEMENT_PROPERTIES[atoms[i].element] || ELEMENT_PROPERTIES.C;
    const z1 = p1.atomicNumber;
    const isMetal1 =
      (z1 >= 21 && z1 <= 30) ||
      (z1 >= 39 && z1 <= 48) ||
      (z1 >= 57 && z1 <= 80);

    for (let j = i + 1; j < atoms.length; j++) {
      const p2 = ELEMENT_PROPERTIES[atoms[j].element] || ELEMENT_PROPERTIES.C;
      const z2 = p2.atomicNumber;
      const isMetal2 =
        (z2 >= 21 && z2 <= 30) ||
        (z2 >= 39 && z2 <= 48) ||
        (z2 >= 57 && z2 <= 80);

      const d = Math.hypot(
        atoms[i].x3d - atoms[j].x3d,
        atoms[i].y3d - atoms[j].y3d,
        atoms[i].z3d - atoms[j].z3d
      );
      if (d < 0.4) continue;

      const rSum = p1.covalentRadius + p2.covalentRadius;
      const tol = isMetal1 || isMetal2 ? 1.28 : 1.22;
      if (d <= rSum * tol) {
        let order: 1 | 2 | 3 = 1;
        // Infer multiple bond order for main-group second-row pairs
        if (!isMetal1 && !isMetal2 && z1 > 1 && z2 > 1) {
          if (d < rSum * 0.80) order = 3;
          else if (d < rSum * 0.90) order = 2;
        }
        bonds.push({
          id: `b_${idx++}`,
          source: atoms[i].id,
          target: atoms[j].id,
          order,
        });
      }
    }
  }
  return bonds;
}

/**
 * Centers 3D atomic coordinates around the geometric centroid and creates a clean 2D projection.
 */
function centerAndProject2D(rawAtoms: Omit<AtomNode, 'x2d' | 'y2d'>[]): AtomNode[] {
  if (rawAtoms.length === 0) return [];
  let cx = 0,
    cy = 0,
    cz = 0;
  for (const a of rawAtoms) {
    cx += a.x3d;
    cy += a.y3d;
    cz += a.z3d;
  }
  cx /= rawAtoms.length;
  cy /= rawAtoms.length;
  cz /= rawAtoms.length;

  return rawAtoms.map((a) => {
    const x3d = a.x3d - cx;
    const y3d = a.y3d - cy;
    const z3d = a.z3d - cz;
    return {
      ...a,
      x3d,
      y3d,
      z3d,
      x2d: x3d * 1.1 + z3d * 0.25,
      y2d: y3d * 1.1 - z3d * 0.2,
    };
  });
}

/**
 * Parses a standard .xyz coordinate file string into AtomNode[] and BondEdge[].
 */
export function parseXYZ(text: string): ParsedXYZResult {
  const lines = text
    .replace(/\r\n/g, '\n')
    .split('\n')
    .map((l) => l.trim());

  const nonEmptyStart = lines.findIndex((l) => l.length > 0);
  if (nonEmptyStart < 0) {
    throw new Error('XYZ file is empty.');
  }

  const countLine = lines[nonEmptyStart];
  const atomCount = parseInt(countLine, 10);
  if (!Number.isFinite(atomCount) || atomCount <= 0 || atomCount > 500) {
    throw new Error(`Invalid XYZ atom count on first line: "${countLine}".`);
  }

  const title = lines[nonEmptyStart + 1] || 'Imported XYZ Structure';
  const rawAtoms: Omit<AtomNode, 'x2d' | 'y2d'>[] = [];

  for (let i = 0; i < atomCount; i++) {
    const lineIdx = nonEmptyStart + 2 + i;
    const line = lines[lineIdx];
    if (!line) {
      throw new Error(`Unexpected end of XYZ file at atom ${i + 1} of ${atomCount}.`);
    }
    const tokens = line.split(/\s+/);
    if (tokens.length < 4) {
      throw new Error(`Malformed coordinate line ${lineIdx + 1}: "${line}".`);
    }
    const el = normalizeElementSymbol(tokens[0]);
    if (!el) {
      throw new Error(`Unsupported or unrecognized element symbol "${tokens[0]}" on line ${lineIdx + 1}.`);
    }
    const x = parseFloat(tokens[1]);
    const y = parseFloat(tokens[2]);
    const z = parseFloat(tokens[3]);
    if (![x, y, z].every(Number.isFinite)) {
      throw new Error(`Invalid numeric coordinates on line ${lineIdx + 1}.`);
    }
    rawAtoms.push({
      id: `${el.toLowerCase()}${i + 1}`,
      element: el,
      x3d: x,
      y3d: y,
      z3d: z,
    });
  }

  const atoms = centerAndProject2D(rawAtoms);
  const bonds = perceiveBondsFrom3D(atoms);
  return { title, atoms, bonds };
}

/**
 * Serializes molecular coordinates to standard .xyz format (Angstroms).
 */
export function exportXYZ(atoms: AtomNode[], title = 'SymmOrbital v3.0 Export'): string {
  const header = `${atoms.length}\n${title.replace(/\r?\n/g, ' ')}`;
  const body = atoms
    .map(
      (a) =>
        `${a.element.padEnd(3, ' ')} ${a.x3d.toFixed(6).padStart(12, ' ')} ${a.y3d
          .toFixed(6)
          .padStart(12, ' ')} ${a.z3d.toFixed(6).padStart(12, ' ')}`
    )
    .join('\n');
  return `${header}\n${body}\n`;
}

/**
 * Parses a Gaussian .cube volumetric file (supports positive Bohr or negative Angstrom axis counts).
 */
export function parseGaussianCube(text: string): GaussianCubeData {
  const lines = text.replace(/\r\n/g, '\n').split('\n');
  if (lines.length < 7) {
    throw new Error('Gaussian .cube file is too short (fewer than 7 header lines).');
  }

  const title = lines[0].trim() || 'Imported Gaussian Cube';
  const comment = lines[1].trim() || '';

  const line3 = lines[2].trim().split(/\s+/).map(Number);
  if (line3.length < 4 || line3.some((v) => !Number.isFinite(v))) {
    throw new Error('Malformed atom count & origin line (line 3) in .cube file.');
  }
  const rawNatoms = Math.trunc(line3[0]);
  const hasDsetHeader = rawNatoms < 0;
  const nAtoms = Math.abs(rawNatoms);
  if (nAtoms === 0 || nAtoms > 500) {
    throw new Error(`Unsupported atom count (${nAtoms}) in .cube file.`);
  }

  const parseAxisLine = (lineStr: string, label: string) => {
    const parts = lineStr.trim().split(/\s+/).map(Number);
    if (parts.length < 4 || parts.some((v) => !Number.isFinite(v))) {
      throw new Error(`Malformed ${label} voxel line in .cube file.`);
    }
    const rawN = Math.trunc(parts[0]);
    const unitScale = rawN < 0 ? 1.0 : BOHR_TO_ANGSTROM;
    return {
      count: Math.abs(rawN),
      vec: [parts[1] * unitScale, parts[2] * unitScale, parts[3] * unitScale] as [number, number, number],
      unitScale,
    };
  };

  const axX = parseAxisLine(lines[3], 'X');
  const axY = parseAxisLine(lines[4], 'Y');
  const axZ = parseAxisLine(lines[5], 'Z');

  if (axX.count < 2 || axY.count < 2 || axZ.count < 2 || axX.count * axY.count * axZ.count > 2_000_000) {
    throw new Error(`Unsupported volumetric grid size (${axX.count}×${axY.count}×${axZ.count}).`);
  }

  const originScale = axX.unitScale;
  const rawOrigin: [number, number, number] = [
    line3[1] * originScale,
    line3[2] * originScale,
    line3[3] * originScale,
  ];

  const rawAtoms: Omit<AtomNode, 'x2d' | 'y2d'>[] = [];
  for (let i = 0; i < nAtoms; i++) {
    const parts = (lines[6 + i] || '').trim().split(/\s+/).map(Number);
    if (parts.length < 5 || parts.some((v) => !Number.isFinite(v))) {
      throw new Error(`Malformed atom line ${7 + i} in .cube file.`);
    }
    const zNum = Math.round(parts[0]);
    const el = ATOMIC_NUMBER_TO_SYMBOL[zNum] || 'C';
    rawAtoms.push({
      id: `${el.toLowerCase()}${i + 1}`,
      element: el,
      x3d: parts[2] * originScale,
      y3d: parts[3] * originScale,
      z3d: parts[4] * originScale,
    });
  }

  // Compute centering shift so both atoms and cube origin stay perfectly registered!
  let cx = 0,
    cy = 0,
    cz = 0;
  for (const a of rawAtoms) {
    cx += a.x3d;
    cy += a.y3d;
    cz += a.z3d;
  }
  cx /= rawAtoms.length;
  cy /= rawAtoms.length;
  cz /= rawAtoms.length;

  const atoms: AtomNode[] = rawAtoms.map((a) => ({
    ...a,
    x3d: a.x3d - cx,
    y3d: a.y3d - cy,
    z3d: a.z3d - cz,
    x2d: (a.x3d - cx) * 1.1 + (a.z3d - cz) * 0.25,
    y2d: (a.y3d - cy) * 1.1 - (a.z3d - cz) * 0.2,
  }));
  const originAngstrom: [number, number, number] = [
    rawOrigin[0] - cx,
    rawOrigin[1] - cy,
    rawOrigin[2] - cz,
  ];
  const bonds = perceiveBondsFrom3D(atoms);

  // Parse volumetric scalar values (x slowest, y middle, z fastest per Gaussian spec)
  const dataStartLine = 6 + nAtoms + (hasDsetHeader ? 1 : 0);
  const totalVoxels = axX.count * axY.count * axZ.count;
  const values = new Float32Array(totalVoxels);
  let cursor = 0;
  let minVal = Infinity;
  let maxVal = -Infinity;

  for (let l = dataStartLine; l < lines.length && cursor < totalVoxels; l++) {
    const trimmed = lines[l].trim();
    if (!trimmed) continue;
    const tokens = trimmed.split(/\s+/);
    for (let t = 0; t < tokens.length && cursor < totalVoxels; t++) {
      const v = parseFloat(tokens[t]);
      const val = Number.isFinite(v) ? v : 0;
      values[cursor++] = val;
      if (val < minVal) minVal = val;
      if (val > maxVal) maxVal = val;
    }
  }

  if (cursor < totalVoxels) {
    throw new Error(
      `Incomplete volumetric cube data: read ${cursor} of ${totalVoxels} grid points.`
    );
  }

  return {
    title,
    comment,
    atoms,
    bonds,
    originAngstrom,
    nx: axX.count,
    ny: axY.count,
    nz: axZ.count,
    vxAngstrom: axX.vec,
    vyAngstrom: axY.vec,
    vzAngstrom: axZ.vec,
    values,
    minVal: Number.isFinite(minVal) ? minVal : 0,
    maxVal: Number.isFinite(maxVal) ? maxVal : 0,
  };
}

/**
 * Exports a uniform 3D volumetric scalar field to a standards-compliant Gaussian .cube file (Bohr units).
 * `values` must be indexed with x slowest, y middle, z fastest: values[(ix * ny + iy) * nz + iz].
 */
export function exportGaussianCube(params: {
  title: string;
  comment: string;
  atoms: AtomNode[];
  originAngstrom: [number, number, number];
  nx: number;
  ny: number;
  nz: number;
  dxAngstrom: number;
  dyAngstrom: number;
  dzAngstrom: number;
  values: Float32Array;
}): string {
  const {
    title,
    comment,
    atoms,
    originAngstrom,
    nx,
    ny,
    nz,
    dxAngstrom,
    dyAngstrom,
    dzAngstrom,
    values,
  } = params;

  const lines: string[] = [];
  lines.push(title.replace(/\r?\n/g, ' '));
  lines.push(comment.replace(/\r?\n/g, ' '));
  lines.push(
    `${String(atoms.length).padStart(5, ' ')} ${(originAngstrom[0] * ANGSTROM_TO_BOHR)
      .toFixed(6)
      .padStart(12, ' ')} ${(originAngstrom[1] * ANGSTROM_TO_BOHR)
      .toFixed(6)
      .padStart(12, ' ')} ${(originAngstrom[2] * ANGSTROM_TO_BOHR)
      .toFixed(6)
      .padStart(12, ' ')}`
  );
  lines.push(
    `${String(nx).padStart(5, ' ')} ${(dxAngstrom * ANGSTROM_TO_BOHR)
      .toFixed(6)
      .padStart(12, ' ')}     0.000000     0.000000`
  );
  lines.push(
    `${String(ny).padStart(5, ' ')}     0.000000 ${(dyAngstrom * ANGSTROM_TO_BOHR)
      .toFixed(6)
      .padStart(12, ' ')}     0.000000`
  );
  lines.push(
    `${String(nz).padStart(5, ' ')}     0.000000     0.000000 ${(dzAngstrom * ANGSTROM_TO_BOHR)
      .toFixed(6)
      .padStart(12, ' ')}`
  );

  for (const a of atoms) {
    const props = ELEMENT_PROPERTIES[a.element] || ELEMENT_PROPERTIES.C;
    const z = props.atomicNumber;
    lines.push(
      `${String(z).padStart(5, ' ')} ${z.toFixed(6).padStart(12, ' ')} ${(
        a.x3d * ANGSTROM_TO_BOHR
      )
        .toFixed(6)
        .padStart(12, ' ')} ${(a.y3d * ANGSTROM_TO_BOHR)
        .toFixed(6)
        .padStart(12, ' ')} ${(a.z3d * ANGSTROM_TO_BOHR).toFixed(6).padStart(12, ' ')}`
    );
  }

  // Write volumetric values (6 values per line, break after each z-pencil per Gaussian spec)
  for (let ix = 0; ix < nx; ix++) {
    for (let iy = 0; iy < ny; iy++) {
      let rowTokens: string[] = [];
      for (let iz = 0; iz < nz; iz++) {
        const val = values[(ix * ny + iy) * nz + iz] ?? 0;
        rowTokens.push(val.toExponential(5).padStart(13, ' '));
        if (rowTokens.length === 6) {
          lines.push(rowTokens.join(' '));
          rowTokens = [];
        }
      }
      if (rowTokens.length > 0) {
        lines.push(rowTokens.join(' '));
      }
    }
  }

  return `${lines.join('\n')}\n`;
}
