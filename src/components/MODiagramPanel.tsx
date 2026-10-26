import React, { useState } from 'react';
import {
  MolecularOrbital,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import {
  formulaToLatex,
  irrepToLatex,
  lcaoExpansionToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { Latex, SvgMathLabel } from './Latex';
import {
  ArrowDown,
  ArrowUp,
  ChevronDown,
  ChevronUp,
  SlidersHorizontal,
  Zap,
} from 'lucide-react';

interface MODiagramPanelProps {
  molecule: MoleculeEntry;
  selectedMoId: string;
  onSelectMo: (mo: MolecularOrbital) => void;
  charge: number;
  onChangeCharge: (charge: number) => void;
  excitedState: boolean;
  onToggleExcitedState: () => void;
  highlightedIrrepId?: string | null;
  onHoverIrrep?: (irrepId: string | null) => void;
  compactWorkbench?: boolean;
}

export const MODiagramPanel: React.FC<MODiagramPanelProps> = ({
  molecule,
  selectedMoId,
  onSelectMo,
  charge,
  onChangeCharge,
  excitedState,
  onToggleExcitedState,
  highlightedIrrepId = null,
  onHoverIrrep,
  compactWorkbench = false,
}) => {
  const [showCorrelations, setShowCorrelations] = useState<boolean>(true);
  const [showHeaderOptions, setShowHeaderOptions] = useState<boolean>(false);
  const [viewMode, setViewMode] = useState<'split' | 'ladder' | 'table'>(
    compactWorkbench ? 'ladder' : 'split'
  );

  const mos = molecule.molecularOrbitals;
  const frags = molecule.fragmentOrbitals;

  const homoMo = mos.find((m) => m.isHOMO) || mos[0];
  const lumoMo = mos.find((m) => m.isLUMO) || mos[mos.length - 1];
  const homoLumoGapEv =
    homoMo && lumoMo && lumoMo.energy > homoMo.energy
      ? Number((lumoMo.energy - homoMo.energy).toFixed(2))
      : 0;
  const photonNm =
    homoLumoGapEv > 0 ? Math.round(1239.84 / homoLumoGapEv) : 0;

  // Group MOs that share the same degenerateGroup or identical energy into visual rows
  interface EnergyRow {
    energy: number;
    orbitals: MolecularOrbital[];
  }

  const rows: EnergyRow[] = [];
  const sortedMos = [...mos].sort((a, b) => b.energy - a.energy); // High energy at top, low energy at bottom
  for (const mo of sortedMos) {
    const existing = rows.find(
      (r) =>
        (mo.degenerateGroup &&
          r.orbitals.some((o) => o.degenerateGroup === mo.degenerateGroup)) ||
        Math.abs(r.energy - mo.energy) < 0.05
    );
    if (existing) {
      existing.orbitals.push(mo);
    } else {
      rows.push({ energy: mo.energy, orbitals: [mo] });
    }
  }

  // Step Up / Down through energy-ordered orbitals
  const energyAscendingMos = [...mos].sort((a, b) => a.energy - b.energy);
  const currentAscIdx = energyAscendingMos.findIndex((m) => m.id === selectedMoId);
  const handleStepOrbital = (delta: 1 | -1) => {
    if (energyAscendingMos.length === 0) return;
    const baseIdx = currentAscIdx >= 0 ? currentAscIdx : 0;
    const nextIdx = Math.max(0, Math.min(energyAscendingMos.length - 1, baseIdx + delta));
    onSelectMo(energyAscendingMos[nextIdx]);
  };

  // SVG Correlation Diagram Geometry
  const SVG_W = 580;
  const SVG_H = 360;
  const PAD_TOP = 32;
  const PAD_BOT = 28;

  const allEnergies = [
    ...mos.map((m) => m.energy),
    ...frags.map((f) => f.energy),
  ];
  const minE = Math.min(...allEnergies) - 2.0;
  const maxE = Math.max(...allEnergies) + 2.0;

  const getMoY = (moId: string): number => {
    const rowIdx = rows.findIndex((r) => r.orbitals.some((o) => o.id === moId));
    if (rowIdx === -1) return SVG_H / 2;
    const rankY =
      PAD_TOP +
      (rowIdx / Math.max(1, rows.length - 1)) * (SVG_H - PAD_TOP - PAD_BOT);
    const mo = mos.find((m) => m.id === moId)!;
    const linearY =
      PAD_TOP +
      ((maxE - mo.energy) / Math.max(1, maxE - minE)) *
        (SVG_H - PAD_TOP - PAD_BOT);
    return Number((0.72 * rankY + 0.28 * linearY).toFixed(1));
  };

  const getFragY = (energy: number): number => {
    const raw =
      PAD_TOP +
      ((maxE - energy) / Math.max(1, maxE - minE)) *
        (SVG_H - PAD_TOP - PAD_BOT);
    return Number(Math.max(PAD_TOP, Math.min(SVG_H - PAD_BOT, raw)).toFixed(1));
  };

  const leftFrags = frags.filter((f) => f.side === 'left');
  const rightFrags = frags.filter((f) => f.side === 'right');
  const activeMo = mos.find((m) => m.id === selectedMoId) || mos[0];

  return (
    <div className="bg-white border border-slate-200 rounded-lg flex flex-col overflow-hidden h-full">
      {/* Header Bar with Collapsible State & Options Toggle (Collapsed by Default) */}
      <div className="px-3.5 py-2.5 bg-slate-50 border-b border-slate-200 flex items-center justify-between gap-2 shrink-0">
        <div className="flex items-center gap-2 min-w-0">
          <h2 className="text-sm font-semibold text-slate-900 truncate">
            MO Energy Diagram
          </h2>
          <span className="text-slate-300" aria-hidden="true">·</span>
          <span className="font-mono text-xs text-blue-700 font-semibold whitespace-nowrap">
            ΔE = {homoLumoGapEv.toFixed(2)} eV
          </span>
        </div>

        <div className="flex items-center gap-1.5 shrink-0">
          <div className="flex items-center gap-1">
            <button
              type="button"
              onClick={() => handleStepOrbital(1)}
              disabled={currentAscIdx >= energyAscendingMos.length - 1}
              title="Step up to next higher-energy MO"
              aria-label="Step up to next higher-energy MO"
              className="p-1 bg-white border border-slate-200 rounded text-slate-700 hover:bg-slate-100 disabled:opacity-40 transition-colors"
            >
              <ArrowUp className="w-3.5 h-3.5" />
            </button>
            <button
              type="button"
              onClick={() => handleStepOrbital(-1)}
              disabled={currentAscIdx <= 0}
              title="Step down to next lower-energy MO"
              aria-label="Step down to next lower-energy MO"
              className="p-1 bg-white border border-slate-200 rounded text-slate-700 hover:bg-slate-100 disabled:opacity-40 transition-colors"
            >
              <ArrowDown className="w-3.5 h-3.5" />
            </button>
          </div>

          <button
            type="button"
            onClick={() => setShowHeaderOptions((v) => !v)}
            aria-expanded={showHeaderOptions}
            className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors flex items-center gap-1 whitespace-nowrap ${
              showHeaderOptions || charge !== 0 || excitedState
                ? 'bg-slate-900 text-white border-slate-900'
                : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-100'
            }`}
          >
            <SlidersHorizontal className="w-3 h-3" />
            <span>State &amp; Options</span>
            {showHeaderOptions ? (
              <ChevronUp className="w-3.5 h-3.5" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5" />
            )}
          </button>
        </div>
      </div>

      {/* Collapsible Electronic State Controls & Frontier Metrics (Collapsed by Default) */}
      {showHeaderOptions && (
        <div className="border-b border-slate-200 bg-slate-50/90 divide-y divide-slate-200/80 shrink-0">
          {/* Row 1: Charge, Excitation, Correlation Lines, and View Mode */}
          <div className="px-3 py-1.5 flex flex-wrap items-center justify-between gap-2">
            <div className="flex flex-wrap items-center gap-1.5">
              {/* Charge State Segmented Control */}
              <div
                className="flex items-center bg-slate-200/70 p-0.5 rounded text-xs"
                role="group"
                aria-label="Molecular Charge State"
              >
                {[
                  { val: 1, label: 'q = +1' },
                  { val: 0, label: 'Neutral' },
                  { val: -1, label: 'q = −1' },
                ].map((item) => (
                  <button
                    key={item.val}
                    type="button"
                    onClick={() => onChangeCharge(item.val)}
                    className={`px-2 py-0.5 text-[11px] font-medium rounded transition-colors whitespace-nowrap shrink-0 ${
                      charge === item.val
                        ? 'bg-white text-slate-900 shadow-xs'
                        : 'text-slate-600 hover:text-slate-900'
                    }`}
                  >
                    {item.label}
                  </button>
                ))}
              </div>

              {/* Photo-Excitation Toggle (HOMO -> LUMO) */}
              <button
                type="button"
                onClick={onToggleExcitedState}
                className={`px-2 py-1 text-[11px] font-medium rounded border transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
                  excitedState
                    ? 'bg-amber-600 text-white border-amber-600'
                    : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-100'
                }`}
              >
                <Zap className="w-3 h-3" />
                <span>
                  {excitedState ? (
                    <>
                      Excited (<Latex math="S_1" />)
                    </>
                  ) : (
                    'HOMO→LUMO'
                  )}
                </span>
              </button>

              <button
                type="button"
                onClick={() => setShowCorrelations((v) => !v)}
                className={`px-2 py-1 text-[11px] font-medium rounded border transition-colors whitespace-nowrap shrink-0 ${
                  showCorrelations
                    ? 'bg-slate-900 text-white border-slate-900'
                    : 'bg-white text-slate-600 border-slate-200 hover:bg-slate-100'
                }`}
              >
                {showCorrelations ? 'Correlations: On' : 'Correlations: Off'}
              </button>
            </div>

            <div className="flex items-center bg-slate-200/70 p-0.5 rounded text-[11px]">
              {(
                [
                  { id: 'ladder', label: 'Diagram' },
                  { id: 'table', label: 'Level List' },
                  { id: 'split', label: 'Split' },
                ] as const
              ).map((m) => (
                <button
                  key={m.id}
                  type="button"
                  onClick={() => setViewMode(m.id)}
                  className={`px-2 py-0.5 font-medium rounded transition-colors ${
                    viewMode === m.id
                      ? 'bg-white text-slate-900 shadow-xs'
                      : 'text-slate-600 hover:text-slate-900'
                  }`}
                >
                  {m.label}
                </button>
              ))}
            </div>
          </div>

          {/* Row 2: Quantitative Frontier Metrics Strip */}
          <div className="px-3 py-1.5 bg-slate-100/60 grid grid-cols-2 sm:grid-cols-4 gap-2 text-xs">
            <div>
              <span className="text-slate-500 block text-[10px]">
                HOMO / SOMO (<Latex math="\varepsilon_{\mathrm{HOMO}}" />)
              </span>
              <span className="font-mono font-semibold text-slate-900 tabular-nums text-[11px]">
                {homoMo ? (
                  <>
                    <Latex math={moLabelToLatex(homoMo.label)} /> ({homoMo.energy.toFixed(2)} eV)
                  </>
                ) : (
                  '—'
                )}
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-[10px]">
                LUMO (<Latex math="\varepsilon_{\mathrm{LUMO}}" />)
              </span>
              <span className="font-mono font-semibold text-slate-900 tabular-nums text-[11px]">
                {lumoMo ? (
                  <>
                    <Latex math={moLabelToLatex(lumoMo.label)} /> ({lumoMo.energy.toFixed(2)} eV)
                  </>
                ) : (
                  '—'
                )}
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-[10px]">
                Frontier Gap (<Latex math="\Delta E" />)
              </span>
              <span className="font-mono font-semibold text-blue-700 tabular-nums text-[11px]">
                {homoLumoGapEv.toFixed(2)} eV {photonNm > 0 ? `(${photonNm} nm)` : ''}
              </span>
            </div>
            <div>
              <span className="text-slate-500 block text-[10px]">
                Valence <Latex math="N_{e^-}" /> · Bond Order
              </span>
              <span className="font-mono font-semibold text-slate-900 tabular-nums text-[11px]">
                {molecule.valenceElectrons} <Latex math="e^-" /> ·{' '}
                <Latex math={`\\mathrm{BO} = ${molecule.bondOrder.toFixed(1)}`} />
              </span>
            </div>
          </div>
        </div>
      )}

      {/* Main Content Area: Interactive SVG Correlation Diagram + Accessible Orbital Selector List */}
      <div
        className={`flex-1 min-h-0 grid grid-cols-1 ${
          viewMode === 'split' ? 'lg:grid-cols-12 divide-y lg:divide-y-0 lg:divide-x' : ''
        } divide-slate-200`}
      >
        {/* Full AO - MO - SALC Correlation Diagram SVG */}
        {(viewMode === 'split' || viewMode === 'ladder') && (
          <div
            className={`${
              viewMode === 'split' ? 'lg:col-span-7' : ''
            } p-2.5 bg-white flex flex-col justify-between min-h-0`}
          >
            <div className="flex items-center justify-between text-[11px] font-medium text-slate-500 px-2 pb-1 border-b border-slate-100 shrink-0">
              <span className="truncate max-w-[34%]">{molecule.leftFragmentTitle}</span>
              <span className="text-slate-900 font-semibold whitespace-nowrap">
                <Latex math={formulaToLatex(molecule.formula)} /> MOs (
                <Latex math={schoenfliesToLatex(molecule.schoenflies)} />)
              </span>
              <span className="truncate text-right max-w-[34%]">{molecule.rightFragmentTitle}</span>
            </div>

            <svg
              viewBox={`0 0 ${SVG_W} ${SVG_H}`}
              className="w-full flex-1 min-h-[260px] select-none block"
              data-testid="mo-correlation-svg"
            >
              {/* Zero Energy Reference Line (Non-bonding threshold ~ -13.6 eV) */}
              <line
                x1={32}
                y1={getFragY(-13.6)}
                x2={SVG_W - 16}
                y2={getFragY(-13.6)}
                stroke="#E2E8F0"
                strokeDasharray="3 3"
                strokeWidth="1"
              />
              <text
                x={36}
                y={getFragY(-13.6) - 4}
                fill="#94A3B8"
                className="font-mono text-[9px]"
              >
                H 1s reference (−13.6 eV)
              </text>

              {/* HOMO-LUMO Gap Shaded Band */}
              {homoMo && lumoMo && (
                <g>
                  <rect
                    x={185}
                    y={getMoY(lumoMo.id)}
                    width={210}
                    height={Math.max(10, getMoY(homoMo.id) - getMoY(lumoMo.id))}
                    fill="#EFF6FF"
                    opacity={0.65}
                    rx={4}
                  />
                  <text
                    x={SVG_W / 2}
                    y={(getMoY(homoMo.id) + getMoY(lumoMo.id)) / 2 + 3}
                    textAnchor="middle"
                    fill="#1D4ED8"
                    className="font-mono text-[10px] font-semibold"
                  >
                    ΔE = {homoLumoGapEv.toFixed(2)} eV
                  </text>
                </g>
              )}

              {/* Correlation Dashed Lines between Fragment Orbitals and MOs */}
              {showCorrelations &&
                rows.map((row) =>
                  row.orbitals.map((mo) => {
                    const moY = getMoY(mo.id);
                    const isSelected = mo.id === selectedMoId;
                    const isIrrepMatch =
                      highlightedIrrepId &&
                      (mo.irrepId === highlightedIrrepId ||
                        mo.irrepSymbol === highlightedIrrepId);
                    return mo.correlatedFragmentIds.map((fragId) => {
                      const frag = frags.find((f) => f.id === fragId);
                      if (!frag) return null;
                      const fragY = getFragY(frag.energy);
                      const startX = frag.side === 'left' ? 110 : SVG_W - 110;
                      const endX = frag.side === 'left' ? 202 : SVG_W - 202;
                      return (
                        <line
                          key={`${mo.id}_${frag.id}`}
                          x1={startX}
                          y1={fragY}
                          x2={endX}
                          y2={moY}
                          stroke={
                            isSelected
                              ? '#2563EB'
                              : isIrrepMatch
                                ? '#0D9488'
                                : '#CBD5E1'
                          }
                          strokeWidth={isSelected || isIrrepMatch ? 1.6 : 1.0}
                          strokeDasharray={isSelected ? 'none' : '3 3'}
                          opacity={isSelected || isIrrepMatch ? 0.9 : 0.65}
                        />
                      );
                    });
                  })
                )}

              {/* Left Fragment Orbitals */}
              {leftFrags.map((frag) => {
                const y = getFragY(frag.energy);
                return (
                  <g key={frag.id} transform={`translate(32, ${y})`}>
                    <line x1={0} y1={0} x2={72} y2={0} stroke="#475569" strokeWidth={2.2} />
                    <SvgMathLabel
                      text={frag.label}
                      x={0}
                      y={-6}
                      textAnchor="start"
                      fill="#334155"
                      className="font-mono text-[10px] font-semibold"
                    />
                    <text x={0} y={12} fill="#64748B" className="font-mono text-[9px]">
                      {frag.energy.toFixed(1)} eV
                    </text>
                  </g>
                );
              })}

              {/* Right Fragment Orbitals (SALCs) */}
              {rightFrags.map((frag) => {
                const y = getFragY(frag.energy);
                return (
                  <g key={frag.id} transform={`translate(${SVG_W - 104}, ${y})`}>
                    <line x1={0} y1={0} x2={72} y2={0} stroke="#0D9488" strokeWidth={2.2} />
                    <SvgMathLabel
                      text={frag.label}
                      x={72}
                      y={-6}
                      textAnchor="end"
                      fill="#0F766E"
                      className="font-mono text-[10px] font-semibold"
                    />
                    <text x={72} y={12} textAnchor="end" fill="#64748B" className="font-mono text-[9px]">
                      {frag.energy.toFixed(1)} eV
                    </text>
                  </g>
                );
              })}

              {/* Center Molecular Orbital Rungs */}
              {rows.map((row) => {
                const y = getMoY(row.orbitals[0].id);
                const count = row.orbitals.length;
                const rungWidth = count === 1 ? 84 : count === 2 ? 62 : 44;
                const gap = count === 1 ? 0 : count === 2 ? 28 : 18;
                const totalW = count * rungWidth + (count - 1) * gap;
                const startX = (SVG_W - totalW) / 2;

                return (
                  <g key={`row_${row.orbitals[0].id}`}>
                    {row.orbitals.map((mo, idx) => {
                      const rx = startX + idx * (rungWidth + gap);
                      const isSelected = mo.id === selectedMoId;
                      const isIrrepHighlighted =
                        highlightedIrrepId &&
                        (mo.irrepId === highlightedIrrepId ||
                          mo.irrepSymbol === highlightedIrrepId);
                      const strokeColor = isSelected
                        ? '#2563EB'
                        : isIrrepHighlighted
                          ? '#0D9488'
                          : mo.isHOMO
                            ? '#0D9488'
                            : mo.isLUMO
                              ? '#D97706'
                              : '#1E293B';

                      return (
                        <g
                          key={mo.id}
                          transform={`translate(${rx}, ${y})`}
                          className="cursor-pointer"
                          onClick={() => onSelectMo(mo)}
                          onMouseEnter={() => onHoverIrrep?.(mo.irrepId)}
                          onMouseLeave={() => onHoverIrrep?.(null)}
                        >
                          {/* Click hit box */}
                          <rect
                            x={-6}
                            y={-14}
                            width={rungWidth + 12}
                            height={28}
                            rx={4}
                            fill={
                              isSelected
                                ? '#DBEAFE'
                                : isIrrepHighlighted
                                  ? '#CCFBF1'
                                  : 'transparent'
                            }
                            stroke={
                              isSelected
                                ? '#2563EB'
                                : isIrrepHighlighted
                                  ? '#0D9488'
                                  : 'transparent'
                            }
                            strokeWidth={1.2}
                          />

                          {/* Energy Rung Line */}
                          <line
                            x1={0}
                            y1={0}
                            x2={rungWidth}
                            y2={0}
                            stroke={strokeColor}
                            strokeWidth={isSelected || isIrrepHighlighted ? 3.5 : 2.5}
                            strokeLinecap="round"
                          />

                          {/* Spin Electrons (Up / Down Arrows) */}
                          {mo.occupancy >= 1 && (
                            <path
                              d={`M ${rungWidth / 2 - 5} 6 L ${rungWidth / 2 - 5} -8 L ${rungWidth / 2 - 8} -4`}
                              fill="none"
                              stroke={isSelected ? '#1D4ED8' : '#0F172A'}
                              strokeWidth={1.8}
                              strokeLinecap="round"
                              strokeLinejoin="round"
                            />
                          )}
                          {mo.occupancy === 2 && (
                            <path
                              d={`M ${rungWidth / 2 + 5} -8 L ${rungWidth / 2 + 5} 6 L ${rungWidth / 2 + 8} 2`}
                              fill="none"
                              stroke={isSelected ? '#1D4ED8' : '#0F172A'}
                              strokeWidth={1.8}
                              strokeLinecap="round"
                              strokeLinejoin="round"
                            />
                          )}

                          {/* Orbital Label */}
                          <SvgMathLabel
                            text={mo.label}
                            x={rungWidth / 2}
                            y={-10}
                            textAnchor="middle"
                            fill={isSelected ? '#1D4ED8' : '#0F172A'}
                            className="font-mono text-[10px] font-bold"
                          />
                          <text
                            x={rungWidth / 2}
                            y={13}
                            textAnchor="middle"
                            fill="#64748B"
                            className="font-mono text-[9px]"
                          >
                            {mo.energy.toFixed(1)} eV
                          </text>
                        </g>
                      );
                    })}
                  </g>
                );
              })}
            </svg>

            {/* Compact LCAO Expansion Strip when in Ladder-Only Mode */}
            {viewMode === 'ladder' && activeMo && (
              <div className="pt-2 border-t border-slate-100 flex items-center justify-between gap-2 text-[11px] bg-slate-50/70 px-2.5 py-1.5 rounded shrink-0 overflow-x-auto">
                <div
                  className="text-slate-900 whitespace-nowrap"
                  data-testid="lcao-latex-equation"
                >
                  <Latex math={lcaoExpansionToLatex(activeMo.label, activeMo.lcao)} />
                </div>
                <button
                  type="button"
                  onClick={() => setViewMode('table')}
                  className="text-[10px] font-medium text-blue-600 hover:underline whitespace-nowrap shrink-0"
                >
                  All Levels →
                </button>
              </div>
            )}
          </div>
        )}

        {/* Interactive Accessible Orbital Level Table */}
        {(viewMode === 'split' || viewMode === 'table') && (
          <div
            className={`${
              viewMode === 'split' ? 'lg:col-span-5' : ''
            } p-3 bg-slate-50/50 flex flex-col justify-between min-h-0`}
          >
            <div className="min-h-0 flex-1 flex flex-col">
              <div className="flex items-center justify-between pb-2 mb-2 border-b border-slate-200 shrink-0">
                <span className="text-xs font-semibold text-slate-800">
                  Orbital Levels (<Latex math="\varepsilon_i" />)
                </span>
                {viewMode === 'table' ? (
                  <button
                    type="button"
                    onClick={() => setViewMode('ladder')}
                    className="text-[11px] font-medium text-blue-600 hover:underline"
                  >
                    ← Back to Diagram
                  </button>
                ) : (
                  <span className="text-[10px] font-mono text-slate-500">
                    Symmetry · Spin · Energy
                  </span>
                )}
              </div>

              <div
                className="space-y-1.5 flex-1 overflow-y-auto pr-1"
                role="listbox"
                aria-label="Molecular Orbital Levels"
              >
                {sortedMos.map((mo) => {
                  const isSelected = mo.id === selectedMoId;
                  const isIrrepMatch =
                    highlightedIrrepId &&
                    (mo.irrepId === highlightedIrrepId ||
                      mo.irrepSymbol === highlightedIrrepId);
                  const spinLatex =
                    mo.occupancy === 2
                      ? '\\uparrow\\downarrow\\;(2e^-)'
                      : mo.occupancy === 1
                        ? '\\uparrow\\;(1e^-)'
                        : '0e^-';
                  const frontierTag = mo.isHOMO
                    ? mo.occupancy === 1
                      ? 'SOMO'
                      : 'HOMO'
                    : mo.isLUMO
                      ? 'LUMO'
                      : null;
                  const interactionLatex =
                    mo.interactionType === 'π*'
                      ? '\\pi^*'
                      : mo.interactionType === 'σ*'
                        ? '\\sigma^*'
                        : mo.interactionType === 'π'
                          ? '\\pi'
                          : mo.interactionType === 'σ'
                            ? '\\sigma'
                            : 'n';

                  return (
                    <button
                      key={mo.id}
                      type="button"
                      onClick={() => onSelectMo(mo)}
                      onMouseEnter={() => onHoverIrrep?.(mo.irrepId)}
                      onMouseLeave={() => onHoverIrrep?.(null)}
                      aria-selected={isSelected}
                      className={`w-full text-left px-2.5 py-1.5 rounded-md border transition-colors flex items-center justify-between gap-2 ${
                        isSelected
                          ? 'bg-blue-50/90 border-blue-600 text-slate-900 shadow-xs'
                          : isIrrepMatch
                            ? 'bg-teal-50/80 border-teal-500 text-slate-900'
                            : 'bg-white border-slate-200 text-slate-700 hover:border-slate-300 hover:bg-slate-50'
                      }`}
                    >
                      <div className="min-w-0">
                        <div className="flex items-center gap-1.5 flex-wrap">
                          <span className="font-mono text-xs font-bold text-slate-900">
                            <Latex math={moLabelToLatex(mo.label)} />
                          </span>
                          <span className="text-slate-300" aria-hidden="true">·</span>
                          <span className="font-mono text-xs font-semibold text-blue-700">
                            <Latex math={irrepToLatex(mo.irrepSymbol)} />
                          </span>
                          {frontierTag && (
                            <>
                              <span className="text-slate-300" aria-hidden="true">·</span>
                              <span
                                className={`font-mono text-[11px] font-bold ${
                                  mo.isHOMO ? 'text-teal-700' : 'text-amber-700'
                                }`}
                              >
                                {frontierTag}
                              </span>
                            </>
                          )}
                        </div>
                        <div className="text-[11px] text-slate-500 flex items-center gap-1.5 mt-0.5 truncate">
                          <span>
                            <Latex math={interactionLatex} /> ({mo.bondingType})
                          </span>
                          <span aria-hidden="true">·</span>
                          <span>
                            {mo.nodalPlanes} node{mo.nodalPlanes === 1 ? '' : 's'}
                          </span>
                          <span aria-hidden="true">·</span>
                          <span className="font-mono">
                            <Latex math={spinLatex} />
                          </span>
                        </div>
                      </div>

                      <div className="text-right shrink-0">
                        <span className="font-mono text-xs font-semibold text-slate-900 tabular-nums block">
                          {mo.energy > 0 ? `+${mo.energy.toFixed(2)}` : mo.energy.toFixed(2)} eV
                        </span>
                      </div>
                    </button>
                  );
                })}
              </div>
            </div>

            {/* Selected Orbital LCAO Expansion Formula in Full KaTeX LaTeX */}
            {activeMo && (
              <div className="mt-2.5 pt-2 border-t border-slate-200 text-xs shrink-0">
                <div className="flex items-center justify-between text-[11px] text-slate-500 mb-1">
                  <span>LCAO Wavefunction Expansion:</span>
                  <span className="font-mono font-semibold text-slate-700">
                    Irrep <Latex math={irrepToLatex(activeMo.irrepSymbol)} /> (
                    <Latex math={schoenfliesToLatex(molecule.schoenflies)} />)
                  </span>
                </div>
                <div
                  className="p-2 bg-white border border-slate-200 rounded text-xs text-slate-900 overflow-x-auto whitespace-nowrap"
                  data-testid="lcao-latex-equation"
                >
                  <Latex math={lcaoExpansionToLatex(activeMo.label, activeMo.lcao)} />
                </div>
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
};
