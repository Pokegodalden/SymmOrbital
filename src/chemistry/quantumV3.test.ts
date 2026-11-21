import test from 'node:test';
import assert from 'node:assert/strict';
import { AtomNode, BondEdge, CANONICAL_MOLECULES } from './moleculesData';
import { POINT_GROUPS } from './symmetryData';
import {
  computeCOOP,
  computeDiatomicSTOOverlap,
  computePDOS,
  runLowdinGeneralizedEHT,
} from './quantumEngineV3';
import {
  computeOpticalSelectionRules,
  detectExactPointGroup3D,
  evaluateElectronicTransition,
  projectMOIrrepsExact,
  wavelengthToVisibleColor,
} from './symmetryV3';
import { runWalshDiagramScan } from './walshEngine';
import {
  exportGaussianCube,
  exportXYZ,
  parseGaussianCube,
  parseXYZ,
} from './fileIO';
import {
  buildVolumetricSurfaceV3,
  compute2DSlicePlane,
  evalElectrostaticPotentialAtPoint,
  evaluateUniformGridForCubeExport,
} from './volumetricV3';

function getPreset(formula: string) {
  const m = CANONICAL_MOLECULES.find((x) => x.formula === formula);
  assert.ok(m, `Preset ${formula} should exist`);
  return m;
}

test('1. Exact diatomic STO overlap integrals match analytical H2 1s-1s benchmark', () => {
  // Exact analytical H(1s)-H(1s) overlap at R = 1.4 Bohr, zeta = 1.3:
  // S = (1 + p + p^2/3) * exp(-p), where p = zeta * R = 1.82
  const R = 1.4;
  const zeta = 1.3;
  const p = zeta * R;
  const exact = (1 + p + (p * p) / 3) * Math.exp(-p);
  const num = computeDiatomicSTOOverlap(1, 0, zeta, 1, 0, zeta, 0, R);
  assert.ok(
    Math.abs(num - exact) < 1e-4,
    `Expected H2 1s-1s overlap ${exact.toFixed(6)}, got ${num.toFixed(6)}`
  );
});

test('2. Full Löwdin generalized EHT satisfies C^T S C = I, H C = S C E, and antibonding destabilization bias', () => {
  const h2Atoms: AtomNode[] = [
    { id: 'h1', element: 'H', x2d: -0.5, y2d: 0, x3d: -0.37, y3d: 0, z3d: 0 },
    { id: 'h2', element: 'H', x2d: 0.5, y2d: 0, x3d: 0.37, y3d: 0, z3d: 0 },
  ];
  const h2Bonds: BondEdge[] = [{ id: 'b1', source: 'h1', target: 'h2', order: 1 }];
  const res = runLowdinGeneralizedEHT(h2Atoms, h2Bonds, 0, false);

  // Verify C^T S C = I
  const n = res.basis.length;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let val = 0;
      for (let a = 0; a < n; a++) {
        for (let b = 0; b < n; b++) {
          val += res.coefficients[a][i] * res.S[a][b] * res.coefficients[b][j];
        }
      }
      const expected = i === j ? 1 : 0;
      assert.ok(Math.abs(val - expected) < 1e-6, `C^T S C [${i}][${j}] = ${val}, expected ${expected}`);
    }
  }

  // Verify H C = S C E
  for (let k = 0; k < n; k++) {
    for (let i = 0; i < n; i++) {
      let hc = 0;
      let sce = 0;
      for (let j = 0; j < n; j++) {
        hc += res.H[i][j] * res.coefficients[j][k];
        sce += res.S[i][j] * res.coefficients[j][k] * res.energies[k];
      }
      assert.ok(Math.abs(hc - sce) < 1e-5, `Secular equation residual too large: ${hc} vs ${sce}`);
    }
  }

  // Verify Pauli repulsion: antibonding destabilization > bonding stabilization relative to H_ii = -13.6 eV
  const stab = -13.6 - res.energies[0];
  const destab = res.energies[1] - -13.6;
  assert.ok(stab > 0 && destab > stab, `Expected destab (${destab.toFixed(2)}) > stab (${stab.toFixed(2)}) > 0`);

  // Verify H2 Mayer bond order ≈ 1.0
  assert.equal(res.mayerBonds.length, 1);
  assert.ok(Math.abs(res.mayerBonds[0].mayerOrder - 1.0) < 0.05, `H2 Mayer BO = ${res.mayerBonds[0].mayerOrder}`);
});

test('3. Mulliken charge sum conservation, water polarity, dipole moment, and Mayer bond order progression', () => {
  const water = getPreset('H2O');
  const resNeutral = runLowdinGeneralizedEHT(water.atoms, water.bonds, 0, false);
  const sumNeutral = resNeutral.mullikenAtoms.reduce((s, a) => s + a.netCharge, 0);
  assert.ok(Math.abs(sumNeutral - 0) < 1e-6, `Neutral H2O charge sum = ${sumNeutral}`);

  // Oxygen must be negative, Hydrogens positive, and dipole > 0
  const oAtom = resNeutral.mullikenAtoms.find((a) => a.element === 'O')!;
  const hAtoms = resNeutral.mullikenAtoms.filter((a) => a.element === 'H');
  assert.ok(oAtom.netCharge < -0.1, `Expected negative O charge, got ${oAtom.netCharge}`);
  assert.ok(hAtoms.every((h) => h.netCharge > 0.05), 'Expected positive H charges');
  assert.ok(resNeutral.dipoleMagnitudeDebye > 0.5, `Expected nonzero H2O dipole, got ${resNeutral.dipoleMagnitudeDebye}`);

  // Cation charge sum = +1
  const resCation = runLowdinGeneralizedEHT(water.atoms, water.bonds, +1, false);
  const sumCation = resCation.mullikenAtoms.reduce((s, a) => s + a.netCharge, 0);
  assert.ok(Math.abs(sumCation - 1) < 1e-6, `H2O+ charge sum = ${sumCation}`);
  assert.equal(resCation.unpairedElectrons, 1, 'H2O+ should be a doublet radical (1 unpaired electron)');

  // Bond order progression: C2H4 (C=C) > H2O (O-H) and N2 (triple) > C2H4 (double)
  const ethene = getPreset('C2H4');
  const n2 = getPreset('N2');
  const resEthene = runLowdinGeneralizedEHT(ethene.atoms, ethene.bonds, 0, false);
  const resN2 = runLowdinGeneralizedEHT(n2.atoms, n2.bonds, 0, false);
  const ccBond = resEthene.mayerBonds.find((b) => b.elementA === 'C' && b.elementB === 'C')!;
  const nnBond = resN2.mayerBonds[0];
  assert.ok(ccBond.mayerOrder > 1.6, `Expected C=C Mayer BO > 1.6, got ${ccBond.mayerOrder}`);
  assert.ok(nnBond.mayerOrder > 2.5, `Expected N≡N Mayer BO > 2.5, got ${nnBond.mayerOrder}`);
  assert.ok(nnBond.mayerOrder > ccBond.mayerOrder, 'Triple bond N2 must exceed double bond C2H4');
});

test('4. PDOS integral and COOP bonding vs antibonding sign separation', () => {
  const n2 = getPreset('N2');
  const eht = runLowdinGeneralizedEHT(n2.atoms, n2.bonds, 0, false);
  const pdos = computePDOS(eht, { sigmaEv: 0.4, numPoints: 200 });
  const dE = (pdos.eMax - pdos.eMin) / (pdos.points.length - 1);
  const area = pdos.points.reduce((s, p) => s + p.total * dE, 0);
  // Area of normalized Gaussians should equal total number of MOs (8 for N2 valence)
  assert.ok(Math.abs(area - eht.basis.length) < 0.25, `PDOS area ${area.toFixed(2)} vs N_MO ${eht.basis.length}`);

  const coop = computeCOOP(eht, n2.bonds, { atomAIndex: 0, atomBIndex: 1 });
  assert.ok(coop.icoopAtFermi > 0.5, `Expected strong positive ICOOP for N2, got ${coop.icoopAtFermi}`);
  assert.ok(coop.moOverlapPopulations[0] > 0, 'Lowest σ bonding MO of N2 must have positive COOP');
  assert.ok(
    coop.moOverlapPopulations[coop.moOverlapPopulations.length - 1] < 0,
    'Highest σ* antibonding MO of N2 must have negative COOP'
  );
});

test('5. Expanded 27-group symmetry library, exact 3D detection, and Wigner irrep projection', () => {
  assert.ok(Object.keys(POINT_GROUPS).length >= 27, `Expected >= 27 point groups, got ${Object.keys(POINT_GROUPS).length}`);
  for (const id of ['C1', 'Ci', 'C2', 'C2h', 'C3h', 'D2', 'D3d', 'D4d', 'D5h', 'S4', 'Th', 'Ih']) {
    assert.ok(POINT_GROUPS[id], `Point group ${id} should be defined`);
  }

  // Test C2h on trans-1,2-difluoroethene
  const transDFE: AtomNode[] = [
    { id: 'c1', element: 'C', x2d: -0.67, y2d: 0, x3d: -0.67, y3d: 0, z3d: 0 },
    { id: 'c2', element: 'C', x2d: 0.67, y2d: 0, x3d: 0.67, y3d: 0, z3d: 0 },
    { id: 'f1', element: 'F', x2d: -1.35, y2d: 1.15, x3d: -1.35, y3d: 1.15, z3d: 0 },
    { id: 'h1', element: 'H', x2d: -1.22, y2d: -0.94, x3d: -1.22, y3d: -0.94, z3d: 0 },
    { id: 'f2', element: 'F', x2d: 1.35, y2d: -1.15, x3d: 1.35, y3d: -1.15, z3d: 0 },
    { id: 'h2', element: 'H', x2d: 1.22, y2d: 0.94, x3d: 1.22, y3d: 0.94, z3d: 0 },
  ];
  const symC2h = detectExactPointGroup3D(transDFE);
  assert.equal(symC2h.pointGroupId, 'C2h');
  assert.equal(symC2h.hasInversion, true);

  // Test D3d on staggered ethane
  const rCH = 1.09;
  const rCC = 1.54;
  const cosTet = Math.cos((19.47 * Math.PI) / 180);
  const sinTet = Math.sin((19.47 * Math.PI) / 180);
  const ethane: AtomNode[] = [
    { id: 'c1', element: 'C', x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: -rCC / 2 },
    { id: 'c2', element: 'C', x2d: 0, y2d: 0, x3d: 0, y3d: 0, z3d: rCC / 2 },
  ];
  for (let k = 0; k < 3; k++) {
    const a1 = (k * 2 * Math.PI) / 3;
    const a2 = a1 + Math.PI / 3;
    ethane.push({
      id: `h1_${k}`,
      element: 'H',
      x2d: 0,
      y2d: 0,
      x3d: rCH * cosTet * Math.cos(a1),
      y3d: rCH * cosTet * Math.sin(a1),
      z3d: -rCC / 2 - rCH * sinTet,
    });
    ethane.push({
      id: `h2_${k}`,
      element: 'H',
      x2d: 0,
      y2d: 0,
      x3d: rCH * cosTet * Math.cos(a2),
      y3d: rCH * cosTet * Math.sin(a2),
      z3d: rCC / 2 + rCH * sinTet,
    });
  }
  const symD3d = detectExactPointGroup3D(ethane);
  assert.equal(symD3d.pointGroupId, 'D3d');

  // Test exact Wigner projection on H2O (C2v)
  const water = getPreset('H2O');
  const symH2O = detectExactPointGroup3D(water.atoms);
  assert.equal(symH2O.pointGroupId, 'C2v');
  const ehtH2O = runLowdinGeneralizedEHT(water.atoms, water.bonds, 0, false);
  const projs = projectMOIrrepsExact(ehtH2O, water.atoms, symH2O);
  assert.equal(projs.length, ehtH2O.energies.length);
  // Lowest orbital (2a1) must be totally symmetric A1
  assert.equal(projs[0].irrepId, 'A1');
});

test('6. Optical selection rules enforce Laporte parity and dipole polarization', () => {
  const water = getPreset('H2O');
  const pgC2v = POINT_GROUPS.C2v;
  const mo1b1 = water.molecularOrbitals.find((m) => m.id === 'mo_1b1')!;
  const mo4a1 = water.molecularOrbitals.find((m) => m.id === 'mo_4a1')!;
  const trans = evaluateElectronicTransition(pgC2v, mo1b1, mo4a1, water.atoms);
  assert.equal(trans.isElectricDipoleAllowed, true);
  assert.ok(trans.allowedPolarizations.includes('x'), `B1 ⊗ A1 = B1 should be x-polarized`);

  // Centrosymmetric D2h ethene: B3u (π) -> B2g (π*) is u -> g Laporte ALLOWED; B3u -> Au would be forbidden
  const ethene = getPreset('C2H4');
  const pgD2h = POINT_GROUPS.D2h;
  const rules = computeOpticalSelectionRules(pgD2h, ethene.molecularOrbitals, ethene.atoms);
  assert.equal(rules.vibrational.isCentrosymmetric, true);
  assert.equal(rules.vibrational.mutualExclusionApplies, true);

  // Visible spectrum wavelength-to-RGB mapping (380–750 nm)
  assert.equal(wavelengthToVisibleColor(250), null, 'UV wavelength should return null');
  assert.equal(wavelengthToVisibleColor(900), null, 'IR wavelength should return null');
  const blue = wavelengthToVisibleColor(455)!;
  assert.ok(blue && blue.bandName === 'Blue' && blue.b > 200 && blue.r === 0);
  const green = wavelengthToVisibleColor(530)!;
  assert.ok(green && green.bandName === 'Green' && green.g > 200);
  const red = wavelengthToVisibleColor(660)!;
  assert.ok(red && red.bandName === 'Red' && red.r > 200 && red.b === 0);
});

test('7. Live Walsh scan on H2O reproduces classic 3a1 bending stabilization and C2v -> Dinfh symmetry change', () => {
  const water = getPreset('H2O');
  const walsh = runWalshDiagramScan(water.atoms, water.bonds, 0, 'bend', 11);
  assert.equal(walsh.coordinate.type, 'bend');
  assert.equal(walsh.frames.length, 11);
  const bentFrame = walsh.frames[0]; // 90 deg
  const midBentFrame = walsh.frames[5]; // 135 deg (EHT bent minimum)
  const nearLinFrame = walsh.frames[walsh.frames.length - 2]; // 171 deg
  const linFrame = walsh.frames[walsh.frames.length - 1]; // 180 deg
  assert.equal(bentFrame.pointGroupId, 'C2v');
  assert.equal(nearLinFrame.pointGroupId, 'C2v');
  assert.equal(linFrame.pointGroupId, 'Dinfh');
  // Classic Walsh rule: 3a1 (MO index 2) stabilizes sharply on bending vs linear π_u, while 1b1 (MO index 3) stays flat
  assert.ok(
    bentFrame.energies[2] < linFrame.energies[2] - 0.5,
    `Expected bent H2O 3a1 (${bentFrame.energies[2].toFixed(2)} eV) < linear πu (${linFrame.energies[2].toFixed(2)} eV)`
  );
  assert.ok(
    Math.abs(bentFrame.energies[3] - linFrame.energies[3]) < 0.05,
    'Non-bonding 1b1 should remain invariant during H2O bending'
  );
  // Bent minimum total electronic energy is strictly lower than linear H2O
  assert.ok(
    midBentFrame.totalElectronicEnergy < linFrame.totalElectronicEnergy,
    `Expected bent H2O E (${midBentFrame.totalElectronicEnergy.toFixed(2)}) < linear E (${linFrame.totalElectronicEnergy.toFixed(2)})`
  );
});

test('8. XYZ and Gaussian .cube round-trip export & import + volumetric density/MEP/2D slice', () => {
  const water = getPreset('H2O');
  const xyzStr = exportXYZ(water.atoms, 'Water Test');
  const parsedXYZ = parseXYZ(xyzStr);
  assert.equal(parsedXYZ.atoms.length, 3);
  assert.equal(parsedXYZ.bonds.length, 2);

  const eht = runLowdinGeneralizedEHT(water.atoms, water.bonds, 0, false);
  const grid = evaluateUniformGridForCubeExport({
    fieldMode: 'density',
    atoms: water.atoms,
    selectedMo: water.molecularOrbitals[3],
    eht,
    resolution: 12,
  });
  const cubeStr = exportGaussianCube({
    title: 'H2O Density',
    comment: 'Unit Test',
    atoms: water.atoms,
    ...grid,
  });
  const parsedCube = parseGaussianCube(cubeStr);
  assert.equal(parsedCube.atoms.length, 3);
  assert.equal(parsedCube.nx, 12);
  assert.equal(parsedCube.values.length, 12 * 12 * 12);
  assert.ok(parsedCube.maxVal > 0.05, 'Density cube should have positive electron density peak');

  // Verify MEP surface & 2D contour slice produce valid geometry and RGBA image
  const mepSurf = buildVolumetricSurfaceV3({
    fieldMode: 'mep',
    atoms: water.atoms,
    eht,
    lobeScale: 1.0,
  });
  assert.ok(mepSurf.hasData, 'MEP surface should produce triangles');
  assert.ok(mepSurf.positive.colors && mepSurf.positive.colors.length === mepSurf.positive.positions.length);

  // H2O 1b1 (O 2pz) has an exact nodal plane on xy (z=0), and strong phase lobes on xz (y=0)
  const nodalSlice = compute2DSlicePlane({
    fieldMode: 'mo',
    atoms: water.atoms,
    selectedMo: water.molecularOrbitals[3],
    eht,
    orientation: 'xy',
    offset: 0,
    showIsolines: true,
    resolution: 32,
  });
  assert.equal(nodalSlice.rgba.length, 32 * 32 * 4);
  assert.ok(nodalSlice.absMax < 1e-4, `Expected exact nodal plane on xy for 1b1, got ${nodalSlice.absMax}`);

  const lobeSlice = compute2DSlicePlane({
    fieldMode: 'mo',
    atoms: water.atoms,
    selectedMo: water.molecularOrbitals[3],
    eht,
    orientation: 'xz',
    offset: 0,
    showIsolines: true,
    resolution: 32,
  });
  assert.equal(lobeSlice.rgba.length, 32 * 32 * 4);
  assert.ok(lobeSlice.absMax > 0.05, `Expected strong 1b1 lobe amplitude on xz slice, got ${lobeSlice.absMax}`);
});

test('9. Quantum DMA MEP reproduces thiophene, benzene, furan, pyrrole, and pyridine electrostatics', () => {
  // 1. Thiophene (C4H4S, experimental microwave C2v geometry in xy plane)
  const thiopheneAtoms: AtomNode[] = [
    { id: 'S1', element: 'S', x2d: 1.46, y2d: 0.0, x3d: 1.46, y3d: 0.0, z3d: 0.0 },
    { id: 'C2', element: 'C', x2d: 0.27, y2d: 1.23, x3d: 0.27, y3d: 1.23, z3d: 0.0 },
    { id: 'C3', element: 'C', x2d: 0.27, y2d: -1.23, x3d: 0.27, y3d: -1.23, z3d: 0.0 },
    { id: 'C4', element: 'C', x2d: -1.0, y2d: 0.71, x3d: -1.0, y3d: 0.71, z3d: 0.0 },
    { id: 'C5', element: 'C', x2d: -1.0, y2d: -0.71, x3d: -1.0, y3d: -0.71, z3d: 0.0 },
    { id: 'H6', element: 'H', x2d: 0.56, y2d: 2.27, x3d: 0.56, y3d: 2.27, z3d: 0.0 },
    { id: 'H7', element: 'H', x2d: 0.56, y2d: -2.27, x3d: 0.56, y3d: -2.27, z3d: 0.0 },
    { id: 'H8', element: 'H', x2d: -1.9, y2d: 1.33, x3d: -1.9, y3d: 1.33, z3d: 0.0 },
    { id: 'H9', element: 'H', x2d: -1.9, y2d: -1.33, x3d: -1.9, y3d: -1.33, z3d: 0.0 },
  ];
  const thiopheneBonds: BondEdge[] = [
    { id: 'b1', source: 'S1', target: 'C2', order: 1 },
    { id: 'b2', source: 'S1', target: 'C3', order: 1 },
    { id: 'b3', source: 'C2', target: 'C4', order: 2 },
    { id: 'b4', source: 'C3', target: 'C5', order: 2 },
    { id: 'b5', source: 'C4', target: 'C5', order: 1 },
    { id: 'b6', source: 'C2', target: 'H6', order: 1 },
    { id: 'b7', source: 'C3', target: 'H7', order: 1 },
    { id: 'b8', source: 'C4', target: 'H8', order: 1 },
    { id: 'b9', source: 'C5', target: 'H9', order: 1 },
  ];
  const ehtThio = runLowdinGeneralizedEHT(thiopheneAtoms, thiopheneBonds, 0, false);
  const vThioPi = evalElectrostaticPotentialAtPoint(ehtThio, thiopheneAtoms, 0.0, 0.0, 1.45);
  const vThioSOut = evalElectrostaticPotentialAtPoint(ehtThio, thiopheneAtoms, 1.46 + 1.8, 0.0, 0.0);
  const vThioSTop = evalElectrostaticPotentialAtPoint(ehtThio, thiopheneAtoms, 1.46, 0.0, 1.75);
  const vThioHAlpha = evalElectrostaticPotentialAtPoint(ehtThio, thiopheneAtoms, 0.85, 3.45, 0.0);
  const vThioHBeta = evalElectrostaticPotentialAtPoint(ehtThio, thiopheneAtoms, -2.9, 2.05, 0.0);

  // Thiophene literature MEP: most negative in ring center, moderately negative on S, positive on H rim
  assert.ok(vThioPi < -0.025, `Thiophene ring center pi should be strongly negative, got ${vThioPi.toFixed(4)}`);
  assert.ok(vThioSOut < -0.005, `Thiophene S outer region should be moderately negative, got ${vThioSOut.toFixed(4)}`);
  assert.ok(vThioPi < vThioSOut, `Thiophene ring center (${vThioPi.toFixed(4)}) must be more negative than S (${vThioSOut.toFixed(4)})`);
  assert.ok(vThioPi < vThioSTop, `Thiophene ring center (${vThioPi.toFixed(4)}) must be more negative than S top (${vThioSTop.toFixed(4)})`);
  assert.ok(vThioHAlpha > 0.015, `Thiophene alpha-H must be positive, got ${vThioHAlpha.toFixed(4)}`);
  assert.ok(vThioHBeta > 0.010, `Thiophene beta-H must be positive, got ${vThioHBeta.toFixed(4)}`);

  const thioSurf = buildVolumetricSurfaceV3({
    fieldMode: 'mep',
    atoms: thiopheneAtoms,
    eht: ehtThio,
    lobeScale: 1.0,
  });
  assert.ok(thioSurf.hasData && thioSurf.mepRange && thioSurf.mepRange.max >= 0.02);

  // 2. Benzene (D6h): strongly negative pi-face, positive H rim
  const rCC = 1.397;
  const rCH = 1.087;
  const benzAtoms: AtomNode[] = [];
  const benzBonds: BondEdge[] = [];
  for (let i = 0; i < 6; i++) {
    const a = (i * Math.PI) / 3;
    benzAtoms.push({ id: `C${i + 1}`, element: 'C', x2d: rCC * Math.cos(a), y2d: rCC * Math.sin(a), x3d: rCC * Math.cos(a), y3d: rCC * Math.sin(a), z3d: 0 });
  }
  for (let i = 0; i < 6; i++) {
    const a = (i * Math.PI) / 3;
    benzAtoms.push({ id: `H${i + 1}`, element: 'H', x2d: (rCC + rCH) * Math.cos(a), y2d: (rCC + rCH) * Math.sin(a), x3d: (rCC + rCH) * Math.cos(a), y3d: (rCC + rCH) * Math.sin(a), z3d: 0 });
    benzBonds.push({ id: `bc${i}`, source: `C${i + 1}`, target: `C${((i + 1) % 6) + 1}`, order: i % 2 === 0 ? 2 : 1 });
    benzBonds.push({ id: `bh${i}`, source: `C${i + 1}`, target: `H${i + 1}`, order: 1 });
  }
  const ehtBenz = runLowdinGeneralizedEHT(benzAtoms, benzBonds, 0, false);
  const vBenzPi = evalElectrostaticPotentialAtPoint(ehtBenz, benzAtoms, 0, 0, 1.6);
  const vBenzH = evalElectrostaticPotentialAtPoint(ehtBenz, benzAtoms, rCC + rCH + 1.15, 0, 0);
  assert.ok(vBenzPi < -0.028, `Benzene pi face should be strongly negative, got ${vBenzPi.toFixed(4)}`);
  assert.ok(vBenzH > 0.012, `Benzene H rim should be positive, got ${vBenzH.toFixed(4)}`);

  // 3. Pyrrole (C4H5N): strongly negative pi-ring, strongly positive N-H
  const pyrroleAtoms: AtomNode[] = [
    { id: 'N1', element: 'N', x2d: 1.12, y2d: 0.0, x3d: 1.12, y3d: 0.0, z3d: 0.0 },
    { id: 'C2', element: 'C', x2d: 0.33, y2d: 1.12, x3d: 0.33, y3d: 1.12, z3d: 0.0 },
    { id: 'C3', element: 'C', x2d: 0.33, y2d: -1.12, x3d: 0.33, y3d: -1.12, z3d: 0.0 },
    { id: 'C4', element: 'C', x2d: -0.98, y2d: 0.71, x3d: -0.98, y3d: 0.71, z3d: 0.0 },
    { id: 'C5', element: 'C', x2d: -0.98, y2d: -0.71, x3d: -0.98, y3d: -0.71, z3d: 0.0 },
    { id: 'H6', element: 'H', x2d: 2.13, y2d: 0.0, x3d: 2.13, y3d: 0.0, z3d: 0.0 },
    { id: 'H7', element: 'H', x2d: 0.69, y2d: 2.14, x3d: 0.69, y3d: 2.14, z3d: 0.0 },
    { id: 'H8', element: 'H', x2d: 0.69, y2d: -2.14, x3d: 0.69, y3d: -2.14, z3d: 0.0 },
    { id: 'H9', element: 'H', x2d: -1.87, y2d: 1.33, x3d: -1.87, y3d: 1.33, z3d: 0.0 },
    { id: 'H10', element: 'H', x2d: -1.87, y2d: -1.33, x3d: -1.87, y3d: -1.33, z3d: 0.0 },
  ];
  const pyrroleBonds: BondEdge[] = [
    { id: 'p1', source: 'N1', target: 'C2', order: 1 },
    { id: 'p2', source: 'N1', target: 'C3', order: 1 },
    { id: 'p3', source: 'C2', target: 'C4', order: 2 },
    { id: 'p4', source: 'C3', target: 'C5', order: 2 },
    { id: 'p5', source: 'C4', target: 'C5', order: 1 },
    { id: 'p6', source: 'N1', target: 'H6', order: 1 },
    { id: 'p7', source: 'C2', target: 'H7', order: 1 },
    { id: 'p8', source: 'C3', target: 'H8', order: 1 },
    { id: 'p9', source: 'C4', target: 'H9', order: 1 },
    { id: 'p10', source: 'C5', target: 'H10', order: 1 },
  ];
  const ehtPyr = runLowdinGeneralizedEHT(pyrroleAtoms, pyrroleBonds, 0, false);
  const vPyrPi = evalElectrostaticPotentialAtPoint(ehtPyr, pyrroleAtoms, 0.0, 0.0, 1.55);
  const vPyrNH = evalElectrostaticPotentialAtPoint(ehtPyr, pyrroleAtoms, 3.25, 0.0, 0.0);
  const vPyrCH = evalElectrostaticPotentialAtPoint(ehtPyr, pyrroleAtoms, -2.85, 2.05, 0.0);
  assert.ok(vPyrPi < -0.04, `Pyrrole pi ring should be deeply negative, got ${vPyrPi.toFixed(4)}`);
  assert.ok(vPyrNH > vPyrCH && vPyrNH > 0.04, `Pyrrole N-H (${vPyrNH.toFixed(4)}) must be strongly positive and > C-H (${vPyrCH.toFixed(4)})`);
});

test('10. Quantum DMA MEP captures N2 quadrupole, chlorobenzene sigma-hole, ethylene pi, and water/cation polarity', () => {
  // 1. N2 (homonuclear diatomic along X): negative lone-pair caps, positive equatorial belt
  const n2 = getPreset('N2');
  const ehtN2 = runLowdinGeneralizedEHT(n2.atoms, n2.bonds, 0, false);
  const vN2Cap = evalElectrostaticPotentialAtPoint(ehtN2, n2.atoms, 2.05, 0, 0);
  const vN2Belt = evalElectrostaticPotentialAtPoint(ehtN2, n2.atoms, 0, 0, 1.6);
  assert.ok(vN2Cap < -0.025, `N2 lone-pair cap should be negative, got ${vN2Cap.toFixed(4)}`);
  assert.ok(vN2Belt > 0.02, `N2 equatorial belt should be positive, got ${vN2Belt.toFixed(4)}`);

  // 2. Ethylene (C2H4): negative pi cloud above C=C, positive H rim
  const c2h4 = getPreset('C2H4');
  const ehtEth = runLowdinGeneralizedEHT(c2h4.atoms, c2h4.bonds, 0, false);
  const vEthPi = evalElectrostaticPotentialAtPoint(ehtEth, c2h4.atoms, 0, 0, 1.6);
  const hEth = c2h4.atoms.find((a) => a.element === 'H')!;
  const vEthH = evalElectrostaticPotentialAtPoint(ehtEth, c2h4.atoms, hEth.x3d * 1.55, hEth.y3d * 1.55, 0);
  assert.ok(vEthPi < -0.02, `C2H4 pi bond should be negative, got ${vEthPi.toFixed(4)}`);
  assert.ok(vEthH > 0.015, `C2H4 H rim should be positive, got ${vEthH.toFixed(4)}`);

  // 3. Chlorobenzene (C6H5Cl): negative aromatic pi-ring, negative Cl equatorial belt, positive Cl sigma-hole
  const rCC = 1.395;
  const rCH = 1.087;
  const rCCl = 1.74;
  const clBenzAtoms: AtomNode[] = [];
  const clBenzBonds: BondEdge[] = [];
  for (let i = 0; i < 6; i++) {
    const a = (i * Math.PI) / 3;
    clBenzAtoms.push({ id: `C${i + 1}`, element: 'C', x2d: rCC * Math.cos(a), y2d: rCC * Math.sin(a), x3d: rCC * Math.cos(a), y3d: rCC * Math.sin(a), z3d: 0 });
  }
  clBenzAtoms.push({ id: 'Cl1', element: 'Cl', x2d: rCC + rCCl, y2d: 0, x3d: rCC + rCCl, y3d: 0, z3d: 0 });
  clBenzBonds.push({ id: 'bcl', source: 'C1', target: 'Cl1', order: 1 });
  for (let i = 0; i < 6; i++) {
    clBenzBonds.push({ id: `bc${i}`, source: `C${i + 1}`, target: `C${((i + 1) % 6) + 1}`, order: i % 2 === 0 ? 2 : 1 });
    if (i > 0) {
      const a = (i * Math.PI) / 3;
      clBenzAtoms.push({ id: `H${i}`, element: 'H', x2d: (rCC + rCH) * Math.cos(a), y2d: (rCC + rCH) * Math.sin(a), x3d: (rCC + rCH) * Math.cos(a), y3d: (rCC + rCH) * Math.sin(a), z3d: 0 });
      clBenzBonds.push({ id: `bh${i}`, source: `C${i + 1}`, target: `H${i}`, order: 1 });
    }
  }
  const ehtClB = runLowdinGeneralizedEHT(clBenzAtoms, clBenzBonds, 0, false);
  const vClBPi = evalElectrostaticPotentialAtPoint(ehtClB, clBenzAtoms, 0, 0, 1.65);
  const vClBelt = evalElectrostaticPotentialAtPoint(ehtClB, clBenzAtoms, rCC + rCCl, 0, 1.75);
  const vClTip = evalElectrostaticPotentialAtPoint(ehtClB, clBenzAtoms, rCC + rCCl + 1.85, 0, 0);
  assert.ok(vClBPi < -0.015, `Chlorobenzene aromatic pi ring should be negative, got ${vClBPi.toFixed(4)}`);
  assert.ok(vClBelt < -0.008, `Chlorobenzene Cl equatorial belt should be negative, got ${vClBelt.toFixed(4)}`);
  assert.ok(vClTip > 0 && vClTip > vClBelt, `Chlorobenzene Cl tip (${vClTip.toFixed(4)}) should exhibit positive sigma-hole vs equatorial belt (${vClBelt.toFixed(4)})`);

  // 4. Water neutral vs H2O+ cation
  const water = getPreset('H2O');
  const ehtH2O = runLowdinGeneralizedEHT(water.atoms, water.bonds, 0, false);
  const ehtH2OCat = runLowdinGeneralizedEHT(water.atoms, water.bonds, +1, false);
  const o = water.atoms[0];
  const h1 = water.atoms[1];
  const h2 = water.atoms[2];
  const mx = 0.5 * (h1.x3d + h2.x3d);
  const my = 0.5 * (h1.y3d + h2.y3d);
  const mz = 0.5 * (h1.z3d + h2.z3d);
  const dl = Math.hypot(o.x3d - mx, o.y3d - my, o.z3d - mz) || 1;
  const lpX = o.x3d + (1.5 * (o.x3d - mx)) / dl;
  const lpY = o.y3d + (1.5 * (o.y3d - my)) / dl;
  const lpZ = o.z3d + (1.5 * (o.z3d - mz)) / dl;
  const vH2OLp = evalElectrostaticPotentialAtPoint(ehtH2O, water.atoms, lpX, lpY, lpZ);
  const vH2OCatLp = evalElectrostaticPotentialAtPoint(ehtH2OCat, water.atoms, lpX, lpY, lpZ);
  assert.ok(vH2OLp < -0.06, `Water O lone pair should be strongly negative, got ${vH2OLp.toFixed(4)}`);
  assert.ok(vH2OCatLp > vH2OLp + 0.15, `H2O+ cation MEP (${vH2OCatLp.toFixed(4)}) must shift strongly positive vs neutral (${vH2OLp.toFixed(4)})`);
});

