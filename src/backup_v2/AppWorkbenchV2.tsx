/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 *
 * ARCHIVED RELEASE (v2.0 Pro 3-Pane Synchronized Workbench)
 * Preserved intact prior to the v3.0 Computational Workbench upgrade.
 * Accessible at any time via the top-right Version Selector Dropdown.
 */

import React, { useEffect, useMemo, useState } from 'react';
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
  PubChemCompoundProfile,
  lookupPubChemCompound,
} from '../chemistry/pubchemService';
import {
  formulaToLatex,
  irrepToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { ChemText, Latex } from '../components/Latex';
import { ChemicalSketcher2D } from '../components/ChemicalSketcher2D';
import { MODiagramPanel } from '../components/MODiagramPanel';
import { OrbitalVisualizer3D } from '../components/OrbitalVisualizer3D';
import { SymmetryGroupPanel } from '../components/SymmetryGroupPanel';
import { DatabaseVerificationPanel } from '../components/DatabaseVerificationPanel';
import { SymmOrbitalLogo } from '../components/SymmOrbitalLogo';
import {
  AppVersion,
  VersionSwitcherDropdown,
} from '../components/VersionSwitcherDropdown';
import {
  ChevronDown,
  ChevronUp,
  Database,
  ExternalLink,
  FlaskConical,
  Moon,
  PanelLeftOpen,
  Scale,
  Search,
  Sun,
  Table2,
} from 'lucide-react';

type PresetCategory = 'all' | 'complexes' | 'hypervalent' | 'diatomics' | 'vsepr' | 'pi' | 'large';
type BottomDockTab = 'symmetry' | 'audit' | 'pubchem' | 'all';

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

interface AppWorkbenchV2Props {
  currentVersion: AppVersion;
  onSelectVersion: (version: AppVersion) => void;
  darkMode: boolean;
  onToggleDarkMode: () => void;
}

export function AppWorkbenchV2({
  currentVersion,
  onSelectVersion,
  darkMode,
  onToggleDarkMode,
}: AppWorkbenchV2Props) {
  const defaultWater = CANONICAL_MOLECULES[0];

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

  const [benchmarksExpanded, setBenchmarksExpanded] = useState<boolean>(true);
  const [sketcherCollapsed, setSketcherCollapsed] = useState<boolean>(false);
  const [presetCategory, setPresetCategory] = useState<PresetCategory>('all');
  const [activeDockTab, setActiveDockTab] = useState<BottomDockTab>('symmetry');

  const [pubchemQuery, setPubchemQuery] = useState<string>('');
  const [pubchemProfile, setPubchemProfile] = useState<PubChemCompoundProfile | null>(null);
  const [pubchemLoading, setPubchemLoading] = useState<boolean>(false);
  const [pubchemError, setPubchemError] = useState<string | null>(null);

  const molecule: MoleculeEntry = useMemo(() => {
    return analyzeMolecularGraph(atoms, bonds, {
      customName: customTitle,
      customCid,
      charge,
      excitedState,
    });
  }, [atoms, bonds, customTitle, customCid, charge, excitedState]);

  const activePointGroup = useMemo(() => {
    const pgKey = pointGroupOverride || molecule.pointGroupId;
    return POINT_GROUPS[pgKey] || POINT_GROUPS.C2v;
  }, [molecule.pointGroupId, pointGroupOverride]);

  const selectedMo: MolecularOrbital = useMemo(() => {
    const found = molecule.molecularOrbitals.find((m) => m.id === selectedMoId);
    if (found) return found;
    const homo = molecule.molecularOrbitals.find((m) => m.isHOMO);
    return homo || molecule.molecularOrbitals[0];
  }, [molecule, selectedMoId]);

  const filteredPresets = useMemo(() => {
    if (presetCategory === 'all') return CANONICAL_MOLECULES;
    const allowed = PRESET_CATEGORIES[presetCategory].formulas;
    return CANONICAL_MOLECULES.filter((m) => allowed.includes(m.formula));
  }, [presetCategory]);

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
    setPubchemError(null);
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
  };

  const handlePubChemLookup = async (queryText: string) => {
    const q = queryText.trim();
    if (!q) return;
    setPubchemLoading(true);
    setPubchemError(null);

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
          // Keep canonical preset active even if live PubChem metadata fetch fails
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
      }
    } catch (err: any) {
      setPubchemError(err?.message || 'Unable to resolve compound from PubChem.');
    } finally {
      setPubchemLoading(false);
    }
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
          <span className="hidden xl:inline text-slate-300" aria-hidden="true">·</span>
          <span className="hidden xl:inline text-xs text-slate-500 truncate max-w-[260px]">
            Pro 3-Pane Workbench (v2.0)
          </span>
        </div>

        {/* Zone 2: Active Molecule & Symmetry Summary */}
        <div className="hidden sm:flex items-center gap-2 text-xs font-mono text-slate-700 min-w-0 truncate">
          <span className="inline-flex items-center gap-1 min-w-0">
            <span>Compound:</span>
            <strong className="text-slate-900 inline-flex items-center gap-1 min-w-0">
              <span className="truncate max-w-[160px] xl:max-w-[220px]">
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
        </div>

        {/* Zone 3: Dark/Light Mode Toggle + Version Dropdown + Reset Action */}
        <div className="flex items-center gap-2 shrink-0">
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
        <section
          aria-label="PubChem Search and Benchmark Molecules"
          className="bg-white border border-slate-200 rounded-lg p-3.5 flex flex-col gap-3 shadow-2xs"
        >
          <div className="flex flex-col md:flex-row items-stretch md:items-center gap-2.5">
            <form
              onSubmit={(e) => {
                e.preventDefault();
                handlePubChemLookup(pubchemQuery);
              }}
              className="flex-1 flex items-center gap-2"
            >
              <label htmlFor="pubchem-hero-omnibar" className="sr-only">
                Search PubChem by Compound Name, Chemical Formula, SMILES, or CID
              </label>
              <div className="relative flex-1">
                <Search className="w-4 h-4 text-slate-400 absolute left-3.5 top-1/2 -translate-y-1/2 pointer-events-none" />
                <input
                  id="pubchem-hero-omnibar"
                  type="text"
                  value={pubchemQuery}
                  onChange={(e) => setPubchemQuery(e.target.value)}
                  placeholder="Search PubChem by compound name, formula, SMILES, or CID (e.g. nitric oxide, CO, benzene, ozone, 962)..."
                  className="w-full pl-10 pr-4 py-2.5 text-sm font-mono bg-slate-50 border border-slate-300 rounded-lg focus:outline-none focus:bg-white focus:border-blue-600 text-slate-900 transition-colors"
                />
              </div>
              <button
                type="submit"
                disabled={pubchemLoading || !pubchemQuery.trim()}
                className="px-5 py-2.5 text-sm font-semibold text-white bg-blue-600 rounded-lg hover:bg-blue-700 disabled:opacity-50 transition-colors flex items-center gap-2 whitespace-nowrap shrink-0"
              >
                <Search className="w-4 h-4" />
                <span>{pubchemLoading ? 'Searching...' : 'Search PubChem'}</span>
              </button>
            </form>

            <button
              type="button"
              onClick={() => setBenchmarksExpanded((prev) => !prev)}
              aria-expanded={benchmarksExpanded}
              className={`px-3.5 py-2.5 text-xs font-semibold rounded-lg border transition-colors flex items-center justify-center gap-2 whitespace-nowrap shrink-0 ${
                benchmarksExpanded
                  ? 'bg-slate-900 text-white border-slate-900'
                  : 'bg-slate-50 text-slate-700 border-slate-300 hover:bg-slate-100'
              }`}
            >
              <FlaskConical className="w-3.5 h-3.5" />
              <span>Benchmark Molecules ({CANONICAL_MOLECULES.length})</span>
              {benchmarksExpanded ? (
                <ChevronUp className="w-4 h-4" />
              ) : (
                <ChevronDown className="w-4 h-4" />
              )}
            </button>
          </div>

          {pubchemError && (
            <div className="px-3 py-2 bg-rose-50 border border-rose-200 rounded text-xs text-rose-800 font-mono">
              ▲ PubChem Query Error: {pubchemError}
            </div>
          )}

          {benchmarksExpanded && (
            <div className="pt-3 border-t border-slate-200 space-y-2.5">
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

                <div className="flex flex-wrap items-center gap-2 text-xs font-mono text-slate-600">
                  <span>
                    Valence: <strong className="text-slate-900">{molecule.valenceElectrons}e⁻</strong>
                  </span>
                  <span className="text-slate-300" aria-hidden="true">·</span>
                  <span>
                    <Latex math={`\\mathrm{BO} = ${molecule.bondOrder.toFixed(1)}`} />
                  </span>
                  <span className="text-slate-300" aria-hidden="true">·</span>
                  <span>
                    <Latex math={`\\mu = ${molecule.dipoleMomentDebye.toFixed(2)}\\;\\mathrm{D}`} />
                  </span>
                  {molecule.pubchemCid > 0 && (
                    <>
                      <span className="text-slate-300" aria-hidden="true">·</span>
                      <a
                        href={`https://pubchem.ncbi.nlm.nih.gov/compound/${molecule.pubchemCid}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="text-blue-600 hover:underline inline-flex items-center gap-1"
                      >
                        <span>PubChem CID {molecule.pubchemCid}</span>
                        <ExternalLink className="w-3 h-3" />
                      </a>
                    </>
                  )}
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
          aria-label="Synchronized 2D Sketcher, MO Energy Diagram, and 3D Orbital Workbench"
          className="grid grid-cols-1 lg:grid-cols-12 gap-4 items-stretch"
        >
          {sketcherCollapsed ? (
            <div className="lg:col-span-1 lg:h-[490px] bg-white border border-slate-200 rounded-lg p-2 flex lg:flex-col items-center justify-between gap-2">
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
            <div className="lg:col-span-4 lg:h-[490px] flex flex-col">
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

          <div
            id="mo-diagram"
            className={`${sketcherCollapsed ? 'lg:col-span-6' : 'lg:col-span-4'} lg:h-[490px] flex flex-col`}
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

          <div
            className={`${sketcherCollapsed ? 'lg:col-span-5' : 'lg:col-span-4'} lg:h-[490px] flex flex-col`}
          >
            <OrbitalVisualizer3D
              molecule={molecule}
              selectedMo={selectedMo}
              pointGroup={activePointGroup}
              highlightedSymmetryElementId={highlightedSymmetryElementId}
            />
          </div>
        </section>

        {/* Contextual Inspector Dock (v2.0 3-Tab Dock) */}
        <section
          id="symmetry-irreps"
          aria-label="Group Theory Character Table and Authoritative Database Verification Dock"
          className="space-y-3"
        >
          <div className="bg-white border border-slate-200 rounded-lg px-4 py-2.5 flex flex-wrap items-center justify-between gap-3">
            <div className="flex items-center gap-2 text-xs">
              <span className="font-semibold text-slate-900">
                Analytical Inspector Dock:
              </span>
              <span className="text-slate-500 hidden sm:inline">
                Hover or click any Character Table row (<Latex math="\Gamma_i" />) to highlight matching orbitals in the MO Diagram above
              </span>
            </div>

            <div
              className="flex flex-wrap items-center gap-1 bg-slate-100 p-0.5 rounded-md"
              role="tablist"
              aria-label="Analytical Inspector Views"
            >
              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'symmetry'}
                onClick={() => setActiveDockTab('symmetry')}
                className={`px-3 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'symmetry'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <Table2 className="w-3.5 h-3.5 text-blue-600" />
                <span>
                  1. Symmetry &amp; Character Table (
                  <Latex math={schoenfliesToLatex(activePointGroup.schoenflies)} />)
                </span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'audit'}
                onClick={() => setActiveDockTab('audit')}
                className={`px-3 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'audit'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <Scale className="w-3.5 h-3.5 text-teal-600" />
                <span>2. NIST WebBook PES &amp; CCCBDB Audit</span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'pubchem'}
                onClick={() => setActiveDockTab('pubchem')}
                className={`px-3 py-1.5 text-xs font-medium rounded transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0 ${
                  activeDockTab === 'pubchem'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                <Database className="w-3.5 h-3.5 text-indigo-600" />
                <span>3. PubChem PUG REST Inspector</span>
              </button>

              <button
                type="button"
                role="tab"
                aria-selected={activeDockTab === 'all'}
                onClick={() => setActiveDockTab('all')}
                className={`px-3 py-1.5 text-xs font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                  activeDockTab === 'all'
                    ? 'bg-white text-slate-900 shadow-xs'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                Show All Expanded
              </button>
            </div>
          </div>

          {(activeDockTab === 'symmetry' || activeDockTab === 'all') && (
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
          )}

          {(activeDockTab === 'audit' || activeDockTab === 'all') && (
            <div id="database-audit">
              <DatabaseVerificationPanel
                molecule={molecule}
                pubchemProfile={pubchemProfile}
                onSelectBenchmarkMolecule={handleSelectPresetMolecule}
              />
            </div>
          )}

          {(activeDockTab === 'pubchem' || activeDockTab === 'all') && (
            <div
              id="pubchem-inspector"
              className="bg-white border border-slate-200 rounded-lg p-4 space-y-3"
            >
              <div className="flex flex-wrap items-center justify-between gap-3">
                <div>
                  <h2 className="text-sm font-semibold text-slate-900">
                    PubChem PUG REST &amp; PUG-View Compound &amp; 3D Conformer Lookup
                  </h2>
                  <p className="text-xs text-slate-500">
                    Fetch live 3D SDF atomic coordinates, IUPAC identifiers, and experimental properties from NCBI PubChem to populate the MO &amp; Symmetry Engine
                  </p>
                </div>
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
                    Active Provenance:{' '}
                    <strong className="text-slate-700">
                      {customTitle
                        ? 'Predicted by Extended Hückel LCAO-MO & Jacobi Eigenvalue Solver'
                        : 'Curated NIST WebBook Photoelectron Spectroscopy & CCCBDB Group Theory Benchmark'}
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
          SymmOrbital v2.0 — Molecular Orbital Theory, LCAO Wavefunction &amp; Point-Group Symmetry Workbench
        </span>
        <span>
          Chemical Data Provenance: NIST Chemistry WebBook PES · NIST CCCBDB · NCBI PubChem PUG REST API
        </span>
      </footer>
    </div>
  );
}
