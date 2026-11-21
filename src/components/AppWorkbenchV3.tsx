/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useMemo, useRef, useState } from 'react';
import {
  AtomNode,
  BondEdge,
  CANONICAL_MOLECULES,
  MolecularOrbital,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import { POINT_GROUPS } from '../chemistry/symmetryData';
import { analyzeMolecularGraph } from '../chemistry/moEngine';
import {
  LowdinEHTResult,
  extractLCAOFromLowdin,
  runLowdinGeneralizedEHT,
} from '../chemistry/quantumEngineV3';
import {
  detectExactPointGroup3D,
  projectMOIrrepsExact,
} from '../chemistry/symmetryV3';
import {
  GaussianCubeData,
  exportGaussianCube,
  exportXYZ,
  parseGaussianCube,
  parseXYZ,
} from '../chemistry/fileIO';
import {
  VolumetricFieldMode,
  evaluateUniformGridForCubeExport,
} from '../chemistry/volumetricV3';
import {
  PubChemCompoundProfile,
  lookupPubChemCompound,
} from '../chemistry/pubchemService';
import {
  formulaToLatex,
  irrepToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { ChemText, Latex } from './Latex';
import { ChemicalSketcher2D } from './ChemicalSketcher2D';
import { MODiagramPanel } from './MODiagramPanel';
import { OrbitalVisualizer3DV3 } from './OrbitalVisualizer3DV3';
import { SymmetryGroupPanel } from './SymmetryGroupPanel';
import { OpticalSelectionRulesPanel } from './OpticalSelectionRulesPanel';
import { PopulationBondingPanel } from './PopulationBondingPanel';
import { WalshDiagramPanel, WalshLiveOverride } from './WalshDiagramPanel';
import { DatabaseVerificationPanel } from './DatabaseVerificationPanel';
import { SymmOrbitalLogo } from './SymmOrbitalLogo';
import { AppVersion, VersionSwitcherDropdown } from './VersionSwitcherDropdown';
import {
  BarChart3,
  ChevronDown,
  ChevronUp,
  Database,
  Download,
  ExternalLink,
  FlaskConical,
  Moon,
  PanelLeftOpen,
  Scale,
  Search,
  Sun,
  Table2,
  TrendingUp,
  Upload,
} from 'lucide-react';

type PresetCategory =
  | 'all'
  | 'complexes'
  | 'hypervalent'
  | 'diatomics'
  | 'vsepr'
  | 'pi'
  | 'large';

type DockTabV3 = 'symmetry' | 'population' | 'walsh' | 'audit' | 'pubchem' | 'all';

const PRESET_CATEGORIES: Record<PresetCategory, { label: string; formulas: string[] }> = {
  all: {
    label: `All Benchmarks (${CANONICAL_MOLECULES.length})`,
    formulas: CANONICAL_MOLECULES.map((m) => m.formula),
  },
  complexes: {
    label: 'd & f Coordination Complexes',
    formulas: ['TiCl4', 'PtCl4', 'CrCl6', 'CeCl6', 'C2H4Cl3Pt'],
  },
  hypervalent: {
    label: 'Hypervalent & Lone-Pair VSEPR',
    formulas: ['XeF6', 'XeF4', 'XeF2', 'XeOF4', 'SF4', 'ClF3', 'PF5', 'SF6'],
  },
  diatomics: {
    label: 'Diatomics & Radicals',
    formulas: ['CO', 'NO', 'O2', 'N2', 'CO2'],
  },
  vsepr: {
    label: 'Main-Group VSEPR',
    formulas: ['H2O', 'NH3', 'CH4', 'BF3', 'PF5', 'SF6', 'XeF2', 'XeOF4'],
  },
  pi: {
    label: 'Conjugated π-Systems',
    formulas: ['C2H4', 'CH2O', 'C6H6', 'C40H56', 'C2H4Cl3Pt'],
  },
  large: {
    label: 'Large molecules',
    formulas: ['C40H56', 'C10H16N5O13P3'],
  },
};

function triggerTextDownload(filename: string, content: string, mime = 'text/plain;charset=utf-8') {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}

interface AppWorkbenchV3Props {
  currentVersion: AppVersion;
  onSelectVersion: (version: AppVersion) => void;
  darkMode: boolean;
  onToggleDarkMode: () => void;
}

export function AppWorkbenchV3({
  currentVersion,
  onSelectVersion,
  darkMode,
  onToggleDarkMode,
}: AppWorkbenchV3Props) {
  const defaultWater = CANONICAL_MOLECULES[0];
  const fileInputRef = useRef<HTMLInputElement | null>(null);

  const [atoms, setAtoms] = useState<AtomNode[]>(defaultWater.atoms);
  const [bonds, setBonds] = useState<BondEdge[]>(defaultWater.bonds);
  const [customTitle, setCustomTitle] = useState<string | undefined>(undefined);
  const [customCid, setCustomCid] = useState<number | undefined>(undefined);
  const [charge, setCharge] = useState<number>(0);
  const [excitedState, setExcitedState] = useState<boolean>(false);
  const [pointGroupOverride, setPointGroupOverride] = useState<string | null>(null);
  const [selectedMoId, setSelectedMoId] = useState<string>('mo_1b1');
  const [highlightedSymmetryElementId, setHighlightedSymmetryElementId] = useState<string | null>(null);
  const [highlightedIrrepId, setHighlightedIrrepId] = useState<string | null>(null);

  // v3.0 Computational & Visual States
  const [fieldMode, setFieldMode] = useState<VolumetricFieldMode>('mo');
  const [eigenvalueSource, setEigenvalueSource] = useState<'curated' | 'lowdin'>('curated');
  const [importedCube, setImportedCube] = useState<GaussianCubeData | null>(null);
  const [fileNotice, setFileNotice] = useState<{ kind: 'ok' | 'err'; text: string } | null>(null);
  const [walshOverride, setWalshOverride] = useState<WalshLiveOverride | null>(null);

  // Workbench UI state
  const [benchmarksExpanded, setBenchmarksExpanded] = useState<boolean>(true);
  const [sketcherCollapsed, setSketcherCollapsed] = useState<boolean>(false);
  const [presetCategory, setPresetCategory] = useState<PresetCategory>('all');
  const [activeDockTab, setActiveDockTab] = useState<DockTabV3>('symmetry');

  // PubChem Lookup state
  const [pubchemQuery, setPubchemQuery] = useState<string>('');
  const [pubchemProfile, setPubchemProfile] = useState<PubChemCompoundProfile | null>(null);
  const [pubchemLoading, setPubchemLoading] = useState<boolean>(false);
  const [pubchemError, setPubchemError] = useState<string | null>(null);

  // 1. Canonical / Graph-analyzed MoleculeEntry (reconstructs 3D coordinates for PubChem 2D-SDF
  //    and custom-sketched molecules while preserving curated NIST benchmarks)
  const baseMolecule: MoleculeEntry = useMemo(() => {
    return analyzeMolecularGraph(atoms, bonds, {
      customName: customTitle,
      customCid,
      charge,
      excitedState,
    });
  }, [atoms, bonds, customTitle, customCid, charge, excitedState]);

  const activeAtoms = walshOverride ? walshOverride.atoms : baseMolecule.atoms;
  const activeBonds = baseMolecule.bonds;

  const isCuratedBenchmark = !customTitle && CANONICAL_MOLECULES.some((m) => m.id === baseMolecule.id);

  // 2. Full Löwdin Generalized EHT Engine (live on reconstructed 3D activeAtoms, including Walsh distortion)
  const lowdinEht: LowdinEHTResult = useMemo(() => {
    return runLowdinGeneralizedEHT(activeAtoms, activeBonds, charge, excitedState);
  }, [activeAtoms, activeBonds, charge, excitedState]);

  // 3. Rigorous 27-Group 3D Symmetry Detection (principal inertia + molecule-aligned elements3D)
  const symV3 = useMemo(() => {
    return detectExactPointGroup3D(activeAtoms);
  }, [activeAtoms]);

  // Resolve active Point Group (user override > Walsh live > curated benchmark or exact 27-group)
  const activePointGroup = useMemo(() => {
    if (pointGroupOverride && POINT_GROUPS[pointGroupOverride]) {
      return POINT_GROUPS[pointGroupOverride];
    }
    if (walshOverride && POINT_GROUPS[walshOverride.frame.pointGroupId]) {
      return POINT_GROUPS[walshOverride.frame.pointGroupId];
    }
    if (isCuratedBenchmark && POINT_GROUPS[baseMolecule.pointGroupId]) {
      return POINT_GROUPS[baseMolecule.pointGroupId];
    }
    return symV3.pointGroup || POINT_GROUPS[baseMolecule.pointGroupId] || POINT_GROUPS.C1;
  }, [pointGroupOverride, walshOverride, isCuratedBenchmark, baseMolecule.pointGroupId, symV3]);

  // 4. Optional Live Löwdin MO Projection (used when user switches to 'lowdin' mode or custom geometry)
  const molecule: MoleculeEntry = useMemo(() => {
    const useLowdinSpectrum = eigenvalueSource === 'lowdin' || Boolean(walshOverride);
    if (!useLowdinSpectrum) {
      return {
        ...baseMolecule,
        dipoleMomentDebye:
          isCuratedBenchmark && charge === 0
            ? baseMolecule.dipoleMomentDebye
            : lowdinEht.dipoleMagnitudeDebye,
      };
    }

    const irrepProjections = projectMOIrrepsExact(lowdinEht, activeAtoms, symV3);
    const irrepCounters: Record<string, number> = {};
    const lowdinMos: MolecularOrbital[] = lowdinEht.energies.map((energy, idx) => {
      const proj = irrepProjections[idx];
      const irrepSym = proj?.irrepSymbol || 'a';
      const irrepId = proj?.irrepId || 'A';
      irrepCounters[irrepSym] = (irrepCounters[irrepSym] || 0) + 1;
      const count = irrepCounters[irrepSym];
      const rawOcc = lowdinEht.occupancies[idx] ?? 0;
      const occ: 0 | 1 | 2 = rawOcc >= 2 ? 2 : rawOcc === 1 ? 1 : 0;
      const isHomo = idx === lowdinEht.homoIndex;
      const isLumo = idx === lowdinEht.lumoIndex;
      const lcao = extractLCAOFromLowdin(lowdinEht, idx);
      const bondType: MolecularOrbital['bondingType'] =
        idx < lowdinEht.homoIndex
          ? 'bonding'
          : isHomo
            ? 'non-bonding'
            : 'antibonding';

      return {
        id: `lowdin_mo_${idx}`,
        label: `${count}${irrepSym.toLowerCase()}`,
        irrepId,
        irrepSymbol: irrepSym.toLowerCase(),
        energy: Number(energy.toFixed(2)),
        occupancy: occ,
        bondingType: bondType,
        interactionType: idx > lowdinEht.homoIndex ? 'σ*' : 'σ',
        nodalPlanes: Math.min(3, Math.max(0, Math.floor(idx / 2))),
        isHOMO: isHomo,
        isLUMO: isLumo,
        correlatedFragmentIds: [],
        lcao,
        description: `Löwdin generalized EHT eigenstate #${idx + 1} (${irrepSym}, ${lcao
          .slice(0, 3)
          .map((t) => `${t.coeff >= 0 ? '+' : ''}${t.coeff.toFixed(2)} ${t.orbital}`)
          .join(' ')}, projection confidence ${Math.round((proj?.confidence ?? 0.5) * 100)}%).`,
      };
    });

    return {
      ...baseMolecule,
      atoms: activeAtoms,
      pointGroupId: activePointGroup.id,
      schoenflies: activePointGroup.schoenflies,
      dipoleMomentDebye: lowdinEht.dipoleMagnitudeDebye,
      molecularOrbitals: lowdinMos.length > 0 ? lowdinMos : baseMolecule.molecularOrbitals,
    };
  }, [
    baseMolecule,
    eigenvalueSource,
    walshOverride,
    isCuratedBenchmark,
    charge,
    lowdinEht,
    activeAtoms,
    symV3,
    activePointGroup,
  ]);

  const selectedMo: MolecularOrbital = useMemo(() => {
    const found = molecule.molecularOrbitals.find((m) => m.id === selectedMoId);
    if (found) return found;
    const homo = molecule.molecularOrbitals.find((m) => m.isHOMO);
    return homo || molecule.molecularOrbitals[0];
  }, [molecule, selectedMoId]);

  const selectedMoIndex = useMemo(() => {
    const idx = molecule.molecularOrbitals.findIndex((m) => m.id === selectedMo?.id);
    return idx >= 0 ? idx : Math.max(0, lowdinEht.homoIndex);
  }, [molecule.molecularOrbitals, selectedMo?.id, lowdinEht.homoIndex]);

  const filteredPresets = useMemo(() => {
    if (presetCategory === 'all') return CANONICAL_MOLECULES;
    const allowed = PRESET_CATEGORIES[presetCategory].formulas;
    return CANONICAL_MOLECULES.filter((m) => allowed.includes(m.formula));
  }, [presetCategory]);

  // Keyboard arrow stepping through MOs
  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent) => {
      if (
        document.activeElement?.tagName === 'INPUT' ||
        document.activeElement?.tagName === 'TEXTAREA' ||
        document.activeElement?.tagName === 'SELECT'
      ) {
        return;
      }
      if (e.key === 'ArrowUp' || e.key === 'ArrowDown') {
        e.preventDefault();
        const asc = [...molecule.molecularOrbitals].sort((a, b) => a.energy - b.energy);
        const curIdx = asc.findIndex((m) => m.id === selectedMo?.id);
        const delta = e.key === 'ArrowUp' ? 1 : -1;
        const nextIdx = Math.max(0, Math.min(asc.length - 1, (curIdx >= 0 ? curIdx : 0) + delta));
        if (asc[nextIdx]) setSelectedMoId(asc[nextIdx].id);
      }
    };
    window.addEventListener('keydown', onKeyDown);
    return () => window.removeEventListener('keydown', onKeyDown);
  }, [molecule.molecularOrbitals, selectedMo?.id]);

  const handleSelectPresetMolecule = (preset: MoleculeEntry) => {
    setAtoms(preset.atoms);
    setBonds(preset.bonds);
    setCustomTitle(undefined);
    setCustomCid(preset.pubchemCid);
    setCharge(0);
    setExcitedState(false);
    setPointGroupOverride(null);
    setHighlightedSymmetryElementId(null);
    setHighlightedIrrepId(null);
    setImportedCube(null);
    setWalshOverride(null);
    setFileNotice(null);
    setPubchemError(null);
    if (fieldMode === 'cube') setFieldMode('mo');
    const homo = preset.molecularOrbitals.find((m) => m.isHOMO) || preset.molecularOrbitals[0];
    if (homo) setSelectedMoId(homo.id);
  };

  const handleGraphChange = (
    nextAtoms: AtomNode[],
    nextBonds: BondEdge[],
    nextTitle?: string
  ) => {
    setAtoms(nextAtoms);
    setBonds(nextBonds);
    setCustomTitle(nextTitle);
    setPointGroupOverride(null);
    setWalshOverride(null);
    setImportedCube(null);
    if (fieldMode === 'cube') setFieldMode('mo');
  };

  const handlePubChemLookup = async (queryText: string) => {
    const q = queryText.trim();
    if (!q) return;
    setPubchemLoading(true);
    setPubchemError(null);
    setFileNotice(null);

    try {
      const lower = q.toLowerCase().replace(/^beta[\s-]+/, 'β-');
      const canonical = CANONICAL_MOLECULES.find(
        (m) =>
          m.name.toLowerCase() === lower ||
          m.name.toLowerCase().replace('β-', 'beta-') === q.toLowerCase() ||
          m.name.toLowerCase().replace('β-', 'beta ') === q.toLowerCase() ||
          (q.toLowerCase() === 'atp' && m.id === 'atp') ||
          (q.toLowerCase() === 'adenosine triphosphate' && m.id === 'atp') ||
          m.iupacName.toLowerCase() === lower ||
          m.formula.toLowerCase() === lower ||
          m.id.toLowerCase() === lower ||
          String(m.pubchemCid) === lower
      );

      if (canonical) {
        handleSelectPresetMolecule(canonical);
        try {
          const profile = await lookupPubChemCompound(
            canonical.pubchemCid ? String(canonical.pubchemCid) : q
          );
          setPubchemProfile(profile);
        } catch {
          // Keep canonical preset active even if live PubChem fetch fails
        }
        return;
      }

      const profile = await lookupPubChemCompound(q);
      setPubchemProfile(profile);

      if (profile.atoms.length > 0) {
        setAtoms(profile.atoms);
        setBonds(profile.bonds);
        setCustomTitle(profile.title || profile.iupacName);
        setCustomCid(profile.cid);
        setCharge(profile.charge || 0);
        setExcitedState(false);
        setPointGroupOverride(null);
        setWalshOverride(null);
        setImportedCube(null);
        if (fieldMode === 'cube') setFieldMode('mo');
      }
    } catch (err: any) {
      setPubchemError(err?.message || 'Unable to resolve compound from PubChem.');
    } finally {
      setPubchemLoading(false);
    }
  };

  // File I/O: Import .xyz or Gaussian .cube
  const handleFileImport = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    e.target.value = '';
    if (!file) return;

    try {
      const text = await file.text();
      const lowerName = file.name.toLowerCase();
      if (lowerName.endsWith('.cube') || lowerName.endsWith('.cub')) {
        const parsed = parseGaussianCube(text);
        setAtoms(parsed.atoms);
        setBonds(parsed.bonds);
        setCustomTitle(parsed.title || file.name.replace(/\.(cube|cub)$/i, ''));
        setCustomCid(undefined);
        setCharge(0);
        setExcitedState(false);
        setPointGroupOverride(null);
        setWalshOverride(null);
        setImportedCube(parsed);
        setFieldMode('cube');
        setFileNotice({
          kind: 'ok',
          text: `Loaded Gaussian .cube "${file.name}" (${parsed.atoms.length} atoms, ${parsed.nx}×${parsed.ny}×${parsed.nz} grid).`,
        });
      } else {
        const parsed = parseXYZ(text);
        setAtoms(parsed.atoms);
        setBonds(parsed.bonds);
        setCustomTitle(parsed.title || file.name.replace(/\.xyz$/i, ''));
        setCustomCid(undefined);
        setCharge(0);
        setExcitedState(false);
        setPointGroupOverride(null);
        setWalshOverride(null);
        setImportedCube(null);
        if (fieldMode === 'cube') setFieldMode('mo');
        setFileNotice({
          kind: 'ok',
          text: `Loaded XYZ "${file.name}" (${parsed.atoms.length} atoms, ${parsed.bonds.length} perceived bonds).`,
        });
      }
    } catch (err: any) {
      setFileNotice({
        kind: 'err',
        text: err?.message || `Failed to parse "${file.name}".`,
      });
    }
  };

  const handleExportXYZ = () => {
    const content = exportXYZ(
      activeAtoms,
      `${molecule.name} (${molecule.formula}) | Point Group ${activePointGroup.schoenflies} | SymmOrbital v3.0`
    );
    triggerTextDownload(`${molecule.formula || 'molecule'}.xyz`, content);
    setFileNotice({
      kind: 'ok',
      text: `Exported ${molecule.formula}.xyz (${activeAtoms.length} atoms, Angstrom).`,
    });
  };

  const handleExportCube = () => {
    const exportField: VolumetricFieldMode = fieldMode === 'cube' ? 'mo' : fieldMode;
    const grid = evaluateUniformGridForCubeExport({
      fieldMode: exportField,
      atoms: activeAtoms,
      selectedMo,
      eht: lowdinEht,
      resolution: 26,
    });
    const content = exportGaussianCube({
      title: `SymmOrbital v3.0 - ${molecule.name} (${exportField.toUpperCase()})`,
      comment: `MO: ${selectedMo?.label || 'N/A'} | Group: ${activePointGroup.schoenflies} | Bohr units`,
      atoms: activeAtoms,
      ...grid,
    });
    triggerTextDownload(`${molecule.formula || 'molecule'}_${exportField}.cube`, content);
    setFileNotice({
      kind: 'ok',
      text: `Exported Gaussian .cube (${exportField.toUpperCase()}, 26³ Bohr grid).`,
    });
  };

  const handleSelectMoByIndex = (idx: number) => {
    const asc = [...molecule.molecularOrbitals].sort((a, b) => a.energy - b.energy);
    const target = asc[Math.max(0, Math.min(asc.length - 1, idx))];
    if (target) setSelectedMoId(target.id);
  };

  return (
    <div className="min-h-screen flex flex-col bg-[#F8FAFC] text-slate-900">
      {/* Strict 3-Zone Top Navigation Bar */}
      <header className="flex flex-wrap lg:flex-nowrap items-center justify-between gap-3 px-5 py-2.5 bg-white border-b border-slate-200 sticky top-0 z-30">
        {/* Zone 1: Brand Identity */}
        <div className="flex items-center gap-2.5 min-w-0 shrink-0">
          <a
            href="#top"
            className="flex items-center gap-2 text-base font-bold tracking-tight text-slate-900 whitespace-nowrap shrink-0"
          >
            <SymmOrbitalLogo className="w-7 h-7 shrink-0" />
            <span>SymmOrbital</span>
          </a>
          <span className="hidden 2xl:inline text-slate-300" aria-hidden="true">·</span>
          <span className="hidden 2xl:inline text-xs text-slate-500 truncate max-w-[230px]">
            Computational MO &amp; Symmetry Suite
          </span>
        </div>

        {/* Zone 2: Live Compound, Group, MO & Dipole Telemetry */}
        <div className="hidden sm:flex items-center gap-2 text-xs font-mono text-slate-700 min-w-0 truncate tabular-nums">
          <span className="inline-flex items-center gap-1 min-w-0">
            <strong className="text-slate-900 inline-flex items-center gap-1 min-w-0">
              <span className="truncate max-w-[140px] xl:max-w-[200px]">
                <ChemText text={molecule.name} />
              </span>
              <span className="shrink-0">
                (<Latex math={formulaToLatex(molecule.formula)} />)
              </span>
            </strong>
          </span>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span>
            Group:{' '}
            <strong className="text-blue-700">
              <Latex math={schoenfliesToLatex(activePointGroup.schoenflies)} />
            </strong>
          </span>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span>
            MO:{' '}
            <strong className="text-teal-700">
              {selectedMo ? (
                <>
                  <Latex math={moLabelToLatex(selectedMo.label)} /> (
                  <Latex math={irrepToLatex(selectedMo.irrepSymbol)} />)
                </>
              ) : (
                '—'
              )}
            </strong>
          </span>
          <span className="hidden xl:inline text-slate-300" aria-hidden="true">·</span>
          <span className="hidden xl:inline">
            <Latex math={`|\\vec{\\mu}| = ${lowdinEht.dipoleMagnitudeDebye.toFixed(2)}\\;\\mathrm{D}`} />
          </span>
        </div>

        {/* Zone 3: Solver Toggle + Theme Toggle + Unified Version Selector Dropdown + Reset */}
        <div className="flex items-center gap-2 shrink-0">
          <div
            className="flex items-center gap-0.5 bg-slate-100 border border-slate-200 p-0.5 rounded-md text-[11px] font-mono shrink-0"
            role="group"
            aria-label="MO Eigenvalue Solver"
          >
            <button
              type="button"
              data-testid="spectrum-curated"
              onClick={() => {
                setEigenvalueSource('curated');
                if (walshOverride) setWalshOverride(null);
              }}
              title={
                isCuratedBenchmark
                  ? 'Use NIST Chemistry WebBook PES calibrated orbital energies'
                  : 'Use calibrated Slater-Koster / ligand-field MO model'
              }
              className={`px-2 py-1 rounded font-medium transition-colors whitespace-nowrap shrink-0 ${
                eigenvalueSource === 'curated' && !walshOverride
                  ? 'bg-white text-slate-900 shadow-2xs'
                  : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              NIST PES
            </button>
            <button
              type="button"
              data-testid="spectrum-lowdin"
              onClick={() => setEigenvalueSource('lowdin')}
              title="Use full generalized Löwdin STO Extended Hückel (S⁻¹/² H S⁻¹/²) eigenvalues"
              className={`px-2 py-1 rounded font-medium transition-colors whitespace-nowrap shrink-0 ${
                eigenvalueSource === 'lowdin' || Boolean(walshOverride)
                  ? 'bg-white text-slate-900 shadow-2xs'
                  : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              Löwdin EHT
            </button>
          </div>

          <button
            type="button"
            onClick={onToggleDarkMode}
            title={darkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'}
            aria-label={darkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'}
            data-testid="theme-toggle-button"
            className="p-1.5 text-xs font-medium text-slate-700 bg-slate-100 border border-slate-200 rounded-md hover:bg-slate-200 transition-colors flex items-center justify-center shrink-0"
          >
            {darkMode ? (
              <Sun className="w-4 h-4 text-amber-400" />
            ) : (
              <Moon className="w-4 h-4 text-slate-600" />
            )}
          </button>

          <VersionSwitcherDropdown
            currentVersion={currentVersion}
            onSelectVersion={onSelectVersion}
          />

          <button
            type="button"
            onClick={() => handleSelectPresetMolecule(defaultWater)}
            className="px-3.5 py-1.5 text-xs font-medium text-white bg-blue-600 rounded-md hover:bg-blue-700 transition-colors whitespace-nowrap shrink-0"
          >
            Reset (<Latex math="\mathrm{H_2O}" />)
          </button>
        </div>
      </header>

      {/* Main Scientific Workbench Container */}
      <main id="top" className="flex-1 max-w-[1440px] w-full mx-auto px-4 sm:px-5 py-3.5 space-y-3.5">
        {/* Top Command Bar: PubChem Search + XYZ/Cube File I/O + Benchmarks Toggle */}
        <section
          aria-label="PubChem Search, File I/O and Benchmark Molecules"
          className="bg-white border border-slate-200 rounded-lg p-3.5 flex flex-col gap-3 shadow-2xs"
        >
          <div className="flex flex-col lg:flex-row items-stretch lg:items-center gap-2">
            {/* PubChem Search Form */}
            <form
              onSubmit={(e) => {
                e.preventDefault();
                handlePubChemLookup(pubchemQuery);
              }}
              className="flex-1 flex items-center gap-2"
            >
              <label htmlFor="pubchem-hero-omnibar-v3" className="sr-only">
                Search PubChem by Compound Name, Chemical Formula, SMILES, or CID
              </label>
              <div className="relative flex-1">
                <Search className="w-4 h-4 text-slate-400 absolute left-3 top-1/2 -translate-y-1/2 pointer-events-none" />
                <input
                  id="pubchem-hero-omnibar-v3"
                  type="text"
                  value={pubchemQuery}
                  onChange={(e) => setPubchemQuery(e.target.value)}
                  placeholder="Search PubChem (e.g. ferrocene, ozone, benzene, Wilkinson's catalyst, 962)..."
                  className="w-full pl-9 pr-3 py-2 text-xs sm:text-sm font-mono bg-slate-50 border border-slate-300 rounded-md focus:outline-none focus:bg-white focus:border-blue-600 text-slate-900 transition-colors"
                />
              </div>
              <button
                type="submit"
                disabled={pubchemLoading || !pubchemQuery.trim()}
                className="px-4 py-2 text-xs font-semibold text-white bg-blue-600 rounded-md hover:bg-blue-700 disabled:opacity-50 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
              >
                <Search className="w-3.5 h-3.5" />
                <span>{pubchemLoading ? 'Searching...' : 'Search'}</span>
              </button>
            </form>

            {/* File I/O Actions (.xyz & Gaussian .cube) + Benchmarks Toggle */}
            <div className="flex flex-wrap items-center gap-1.5 shrink-0">
              <input
                ref={fileInputRef}
                type="file"
                accept=".xyz,.cube,.cub,text/plain"
                onChange={handleFileImport}
                className="hidden"
                data-testid="file-import-input"
              />
              <button
                type="button"
                onClick={() => fileInputRef.current?.click()}
                data-testid="import-file-button"
                title="Import Cartesian .xyz or Gaussian .cube volumetric grid"
                className="px-2.5 py-2 text-xs font-medium text-slate-700 bg-slate-50 border border-slate-300 rounded-md hover:bg-slate-100 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
              >
                <Upload className="w-3.5 h-3.5 text-blue-600" />
                <span>Import .xyz/.cube</span>
              </button>

              <button
                type="button"
                onClick={handleExportXYZ}
                data-testid="export-xyz-button"
                title="Export active 3D coordinates as .xyz"
                className="px-2.5 py-2 text-xs font-medium text-slate-700 bg-slate-50 border border-slate-300 rounded-md hover:bg-slate-100 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
              >
                <Download className="w-3.5 h-3.5 text-slate-600" />
                <span>.xyz</span>
              </button>

              <button
                type="button"
                onClick={handleExportCube}
                data-testid="export-cube-button"
                title="Export active 3D volumetric field as Gaussian .cube (Bohr units)"
                className="px-2.5 py-2 text-xs font-medium text-slate-700 bg-slate-50 border border-slate-300 rounded-md hover:bg-slate-100 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
              >
                <Download className="w-3.5 h-3.5 text-teal-600" />
                <span>.cube</span>
              </button>

              <button
                type="button"
                onClick={() => setBenchmarksExpanded((prev) => !prev)}
                aria-expanded={benchmarksExpanded}
                className={`px-3 py-2 text-xs font-semibold rounded-md border transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  benchmarksExpanded
                    ? 'bg-slate-900 text-white border-slate-900'
                    : 'bg-slate-50 text-slate-700 border-slate-300 hover:bg-slate-100'
                }`}
              >
                <FlaskConical className="w-3.5 h-3.5" />
                <span>Benchmarks ({CANONICAL_MOLECULES.length})</span>
                {benchmarksExpanded ? (
                  <ChevronUp className="w-3.5 h-3.5" />
                ) : (
                  <ChevronDown className="w-3.5 h-3.5" />
                )}
              </button>
            </div>
          </div>

          {pubchemError && (
            <div className="px-3 py-2 bg-rose-50 border border-rose-200 rounded text-xs text-rose-800 font-mono">
              ▲ PubChem Query Error: {pubchemError}
            </div>
          )}

          {fileNotice && (
            <div
              className={`px-3 py-1.5 rounded text-xs font-mono border flex items-center justify-between gap-2 ${
                fileNotice.kind === 'ok'
                  ? 'bg-teal-50 border-teal-200 text-teal-800'
                  : 'bg-rose-50 border-rose-200 text-rose-800'
              }`}
            >
              <span>{fileNotice.text}</span>
              <button
                type="button"
                onClick={() => setFileNotice(null)}
                className="text-[11px] underline shrink-0"
              >
                Dismiss
              </button>
            </div>
          )}

          {/* Collapsible Benchmark Molecules Drawer */}
          {benchmarksExpanded && (
            <div className="pt-2.5 border-t border-slate-200 space-y-2.5">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div
                  className="flex flex-wrap items-center gap-1 bg-slate-100 p-0.5 rounded-md"
                  role="tablist"
                  aria-label="Molecule Categories"
                >
                  {(Object.keys(PRESET_CATEGORIES) as PresetCategory[]).map((catKey) => {
                    const isCat = presetCategory === catKey;
                    return (
                      <button
                        key={catKey}
                        type="button"
                        role="tab"
                        aria-selected={isCat}
                        onClick={() => setPresetCategory(catKey)}
                        className={`px-2.5 py-1 text-xs font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                          isCat
                            ? 'bg-white text-slate-900 shadow-xs'
                            : 'text-slate-600 hover:text-slate-900'
                        }`}
                      >
                        <ChemText text={PRESET_CATEGORIES[catKey].label} />
                      </button>
                    );
                  })}
                </div>

                {/* Valence & Bond Order Readout */}
                <div className="flex flex-wrap items-center gap-2 text-xs font-mono text-slate-600">
                  <span className="tabular-nums">
                    Valence: <strong className="text-slate-900">{molecule.valenceElectrons}e⁻</strong>
                  </span>
                  <span className="text-slate-300" aria-hidden="true">·</span>
                  <span className="tabular-nums">
                    <Latex math={`\\mathrm{BO} = ${molecule.bondOrder.toFixed(1)}`} />
                  </span>
                </div>
              </div>

              <div
                className="flex flex-wrap items-center gap-1.5"
                role="group"
                aria-label="Benchmark Molecule Presets"
              >
                {filteredPresets.map((preset) => {
                  const isCurrent = molecule.formula === preset.formula && !customTitle;
                  return (
                    <button
                      key={preset.id}
                      type="button"
                      onClick={() => handleSelectPresetMolecule(preset)}
                      className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors whitespace-nowrap shrink-0 ${
                        isCurrent
                          ? 'bg-slate-900 text-white border-slate-900 shadow-xs'
                          : 'bg-slate-50 text-slate-700 border-slate-200 hover:bg-slate-100 hover:text-slate-900'
                      }`}
                    >
                      <Latex math={formulaToLatex(preset.formula)} /> (<ChemText text={preset.name} />) ·{' '}
                      <Latex math={schoenfliesToLatex(preset.schoenflies)} />
                    </button>
                  );
                })}
              </div>
            </div>
          )}
        </section>

        {/* Primary Synchronized Above-the-Fold 3-Pane Scientific Stage */}
        <section
          id="sketcher-3d"
          aria-label="Synchronized 2D Sketcher, MO Energy Diagram, and 3D Volumetric Field Workbench"
          className="grid grid-cols-1 lg:grid-cols-12 gap-4 items-stretch"
        >
          {/* Column 1: 2D Chemical Sketcher */}
          {sketcherCollapsed ? (
            <div className="lg:col-span-1 lg:h-[510px] bg-white border border-slate-200 rounded-lg p-2 flex lg:flex-col items-center justify-between gap-2">
              <button
                type="button"
                onClick={() => setSketcherCollapsed(false)}
                className="px-2.5 py-2 text-xs font-medium text-blue-700 bg-blue-50 border border-blue-200 rounded hover:bg-blue-100 transition-colors flex items-center gap-1.5"
                title="Expand 2D Chemical Structure Sketcher"
              >
                <PanelLeftOpen className="w-4 h-4 shrink-0" />
                <span className="lg:hidden">Expand 2D Sketcher</span>
              </button>

              <div className="text-xs font-mono text-slate-600 text-center">
                <strong className="text-slate-900 block">2D Sketcher</strong>
                <span>
                  <Latex math={formulaToLatex(molecule.formula)} /> ({atoms.length} atoms)
                </span>
              </div>

              <button
                type="button"
                onClick={() => setSketcherCollapsed(false)}
                className="px-2 py-1 text-[11px] font-medium text-slate-600 hover:text-slate-900 underline"
              >
                Edit 2D
              </button>
            </div>
          ) : (
            <div className="lg:col-span-4 lg:h-[510px] flex flex-col">
              <ChemicalSketcher2D
                atoms={atoms}
                bonds={bonds}
                onChangeGraph={handleGraphChange}
                onTriggerPubChemSearch={handlePubChemLookup}
                isSearchingPubChem={pubchemLoading}
                onCollapseSketcher={() => setSketcherCollapsed(true)}
                hideBottomSearch
              />
            </div>
          )}

          {/* Column 2: Interactive MO Energy Diagram */}
          <div
            id="mo-diagram"
            className={`${sketcherCollapsed ? 'lg:col-span-6' : 'lg:col-span-4'} lg:h-[510px] flex flex-col`}
          >
            <MODiagramPanel
              molecule={molecule}
              selectedMoId={selectedMo?.id || ''}
              onSelectMo={(mo) => setSelectedMoId(mo.id)}
              charge={charge}
              onChangeCharge={setCharge}
              excitedState={excitedState}
              onToggleExcitedState={() => setExcitedState((prev) => !prev)}
              highlightedIrrepId={highlightedIrrepId}
              onHoverIrrep={setHighlightedIrrepId}
              compactWorkbench={!sketcherCollapsed}
            />
          </div>

          {/* Column 3: Upgraded v3.0 3D Volumetric & 2D Contour Slice Visualizer */}
          <div
            className={`${sketcherCollapsed ? 'lg:col-span-5' : 'lg:col-span-4'} lg:h-[510px] flex flex-col`}
          >
            <OrbitalVisualizer3DV3
              molecule={molecule}
              displayAtoms={activeAtoms}
              selectedMo={selectedMo}
              pointGroup={activePointGroup}
              alignedSymmetryElements={
                walshOverride || !isCuratedBenchmark ? symV3.alignedElements3D : undefined
              }
              eht={lowdinEht}
              fieldMode={fieldMode}
              onChangeFieldMode={setFieldMode}
              importedCube={importedCube}
              highlightedSymmetryElementId={highlightedSymmetryElementId}
              walshBanner={
                walshOverride
                  ? `${walshOverride.coordLabel} = ${walshOverride.valueFormatted} (${walshOverride.frame.schoenflies})`
                  : null
              }
              onInspectOpenShellO2={() => {
                const o2 = CANONICAL_MOLECULES.find((m) => m.formula === 'O2');
                if (o2) {
                  handleSelectPresetMolecule(o2);
                  setFieldMode('spin');
                }
              }}
              onIonizeRadical={() => {
                setCharge((c) => c + 1);
                setFieldMode('spin');
              }}
            />
          </div>
        </section>

        {/* v3.0 5-Tab Analytical Inspector Dock */}
        <section
          id="symmetry-irreps"
          aria-label="v3.0 Analytical Computational Chemistry Inspector Dock"
          className="space-y-3"
        >
          <div className="bg-white border border-slate-200 rounded-lg px-4 py-2.5 flex flex-wrap items-center justify-between gap-2.5">
            <span className="text-xs font-semibold text-slate-900 whitespace-nowrap">
              Analytical Inspector:
            </span>

            <div
              className="flex flex-wrap items-center gap-1 bg-slate-100 p-0.5 rounded-md"
              role="tablist"
              aria-label="v3.0 Analytical Views"
            >
              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'symmetry'}
                data-testid="dock-tab-symmetry"
                onClick={() => setActiveDockTab('symmetry')}
                className={`px-2.5 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'symmetry'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <Table2 className="w-3.5 h-3.5 text-blue-600" />
                <span>
                  1. Symmetry &amp; Optical Rules (
                  <Latex math={schoenfliesToLatex(activePointGroup.schoenflies)} />)
                </span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'population'}
                data-testid="dock-tab-population"
                onClick={() => setActiveDockTab('population')}
                className={`px-2.5 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'population'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <BarChart3 className="w-3.5 h-3.5 text-teal-600" />
                <span>2. Mulliken · Mayer · PDOS · COOP</span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'walsh'}
                data-testid="dock-tab-walsh"
                onClick={() => setActiveDockTab('walsh')}
                className={`px-2.5 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'walsh'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <TrendingUp className="w-3.5 h-3.5 text-blue-600" />
                <span>3. Walsh Diagram</span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'audit'}
                data-testid="dock-tab-audit"
                onClick={() => setActiveDockTab('audit')}
                className={`px-2.5 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'audit'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <Scale className="w-3.5 h-3.5 text-teal-600" />
                <span>4. NIST PES Audit</span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'pubchem'}
                data-testid="dock-tab-pubchem"
                onClick={() => setActiveDockTab('pubchem')}
                className={`px-2.5 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'pubchem'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <Database className="w-3.5 h-3.5 text-indigo-600" />
                <span>5. PubChem</span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'all'}
                onClick={() => setActiveDockTab('all')}
                className={`px-2.5 py-1.5 text-xs font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeDockTab === 'all'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                All
              </button>
            </div>
          </div>

          {/* Tab 1: 27-Group Symmetry Character Table + Optical Selection Rules */}
          {(activeDockTab === 'symmetry' || activeDockTab === 'all') && (
            <div className="space-y-3.5">
              <SymmetryGroupPanel
                molecule={molecule}
                pointGroup={activePointGroup}
                selectedMo={selectedMo}
                highlightedSymmetryElementId={highlightedSymmetryElementId}
                onHighlightSymmetryElement={setHighlightedSymmetryElementId}
                onSelectPointGroupOverride={(pgId) => setPointGroupOverride(pgId)}
                highlightedIrrepId={highlightedIrrepId}
                onHoverIrrep={setHighlightedIrrepId}
                onSelectMo={(mo) => setSelectedMoId(mo.id)}
              />

              <OpticalSelectionRulesPanel
                pointGroup={activePointGroup}
                molecularOrbitals={molecule.molecularOrbitals}
                atoms={activeAtoms}
                selectedMo={selectedMo}
                onSelectMo={(mo) => setSelectedMoId(mo.id)}
              />
            </div>
          )}

          {/* Tab 2: Mulliken Charges, Mayer Bond Orders, PDOS & COOP */}
          {(activeDockTab === 'population' || activeDockTab === 'all') && (
            <PopulationBondingPanel
              eht={lowdinEht}
              atoms={activeAtoms}
              bonds={activeBonds}
              molecularOrbitals={molecule.molecularOrbitals}
              selectedMoId={selectedMo?.id || ''}
              onSelectMoIndex={handleSelectMoByIndex}
              darkMode={darkMode}
            />
          )}

          {/* Tab 3: Interactive Walsh Diagram */}
          {(activeDockTab === 'walsh' || activeDockTab === 'all') && (
            <WalshDiagramPanel
              atoms={baseMolecule.atoms}
              bonds={activeBonds}
              charge={charge}
              moleculeName={molecule.name}
              moleculeFormula={molecule.formula}
              selectedMoIndex={selectedMoIndex}
              onWalshGeometryOverride={setWalshOverride}
              onSelectMoIndex={handleSelectMoByIndex}
              darkMode={darkMode}
            />
          )}

          {/* Tab 4: NIST WebBook PES & CCCBDB Audit */}
          {(activeDockTab === 'audit' || activeDockTab === 'all') && (
            <div id="database-audit">
              <DatabaseVerificationPanel
                molecule={molecule}
                pubchemProfile={pubchemProfile}
                onSelectBenchmarkMolecule={handleSelectPresetMolecule}
              />
            </div>
          )}

          {/* Tab 5: PubChem PUG REST Inspector */}
          {(activeDockTab === 'pubchem' || activeDockTab === 'all') && (
            <div
              id="pubchem-inspector"
              className="bg-white border border-slate-200 rounded-lg p-4 space-y-3"
            >
              <div>
                <h2 className="text-sm font-semibold text-slate-900">
                  PubChem PUG REST Compound &amp; 3D Conformer Inspector
                </h2>
                <p className="text-xs text-slate-500">
                  Live 3D SDF atomic coordinates, IUPAC identifiers, and physicochemical telemetry fetched from NCBI PubChem
                </p>
              </div>

              {pubchemProfile ? (
                <div className="grid grid-cols-1 md:grid-cols-4 gap-4 pt-2 border-t border-slate-100 text-xs">
                  <div>
                    <span className="text-slate-500 block text-[11px]">Compound / IUPAC Name</span>
                    <strong className="text-slate-900 font-medium block mt-0.5">
                      {pubchemProfile.title} ({pubchemProfile.iupacName})
                    </strong>
                    <span className="text-[11px] text-slate-500 font-mono mt-1 block">
                      3D Geometry:{' '}
                      {pubchemProfile.conformerType === '3D-SDF'
                        ? 'PubChem MMFF94s 3D-SDF Conformer'
                        : '2D-SDF → 3D VSEPR Reconstructed'}
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-500 block text-[11px]">Formula &amp; Molecular Weight</span>
                    <span className="font-mono font-semibold text-slate-900 tabular-nums block mt-0.5">
                      <Latex math={formulaToLatex(pubchemProfile.molecularFormula)} /> ·{' '}
                      {pubchemProfile.molecularWeight.toFixed(2)} g/mol
                    </span>
                    <span className="font-mono text-[11px] text-slate-500 block mt-0.5 truncate">
                      InChIKey: {pubchemProfile.inchikey}
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-500 block text-[11px]">Canonical SMILES &amp; Polarity</span>
                    <span className="font-mono font-semibold text-blue-700 block mt-0.5 truncate">
                      {pubchemProfile.canonicalSmiles}
                    </span>
                    <span className="font-mono text-[11px] text-slate-600 block mt-0.5 tabular-nums">
                      TPSA: {pubchemProfile.tpsa ?? 0} <Latex math="\text{\AA}^2" /> · XLogP:{' '}
                      {pubchemProfile.xlogp ?? 'N/A'}
                    </span>
                  </div>
                  <div className="flex flex-col justify-between items-start md:items-end">
                    <span className="text-slate-500 text-[11px]">Authoritative Record</span>
                    <a
                      href={pubchemProfile.pubchemRecordUrl}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="px-3 py-1.5 bg-blue-50 text-blue-700 border border-blue-200 rounded font-mono text-xs font-medium hover:bg-blue-100 transition-colors inline-flex items-center gap-1.5 whitespace-nowrap"
                    >
                      <span>View PubChem CID {pubchemProfile.cid}</span>
                      <ExternalLink className="w-3.5 h-3.5" />
                    </a>
                  </div>
                </div>
              ) : (
                <div className="pt-2 border-t border-slate-100 flex flex-wrap items-center justify-between gap-2 text-xs text-slate-500">
                  <span>
                    Active Engine:{' '}
                    <strong className="text-slate-700">
                      Löwdin Generalized EHT (<Latex math="\mathbf{S}^{-1/2}\mathbf{H}\mathbf{S}^{-1/2}" />) + 27-Group 3D Symmetry Projection
                    </strong>
                  </span>
                  {molecule.pubchemCid > 0 && (
                    <a
                      href={`https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.pubchemCid}`}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="font-mono text-blue-600 hover:underline inline-flex items-center gap-1"
                    >
                      <span>PubChem Reference Record CID {molecule.pubchemCid}</span>
                      <ExternalLink className="w-3 h-3" />
                    </a>
                  )}
                </div>
              )}
            </div>
          )}
        </section>
      </main>

      <footer className="px-5 py-3.5 bg-white border-t border-slate-200 text-xs text-slate-500 flex flex-wrap items-center justify-between gap-4">
        <span>
          SymmOrbital v3.0 — Generalized Löwdin EHT, 27-Group Symmetry, Volumetric MEP/Slice &amp; Walsh Correlation Workbench
        </span>
        <span>
          Data Provenance: NIST Chemistry WebBook PES · NIST CCCBDB · Clementi–Raimondi STOs · NCBI PubChem PUG REST
        </span>
      </footer>
    </div>
  );
}
