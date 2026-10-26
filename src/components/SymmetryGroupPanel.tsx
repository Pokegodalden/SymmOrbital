import React, { useEffect, useState } from 'react';
import {
  POINT_GROUPS,
  PointGroupData,
  decomposeReducibleRepresentation,
} from '../chemistry/symmetryData';
import {
  MolecularOrbital,
  MoleculeEntry,
} from '../chemistry/moleculesData';
import {
  basisFunctionToLatex,
  classSymbolToLatex,
  irrepToLatex,
  moLabelToLatex,
  reducibleLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { Latex } from './Latex';
import { ChevronDown, ChevronUp, SlidersHorizontal } from 'lucide-react';

interface SymmetryGroupPanelProps {
  molecule: MoleculeEntry;
  pointGroup: PointGroupData;
  selectedMo: MolecularOrbital;
  highlightedSymmetryElementId: string | null;
  onHighlightSymmetryElement: (id: string | null) => void;
  onSelectPointGroupOverride?: (pgId: string) => void;
  highlightedIrrepId?: string | null;
  onHoverIrrep?: (irrepId: string | null) => void;
  onSelectMo?: (mo: MolecularOrbital) => void;
}

export const SymmetryGroupPanel: React.FC<SymmetryGroupPanelProps> = ({
  molecule,
  pointGroup,
  selectedMo,
  highlightedSymmetryElementId,
  onHighlightSymmetryElement,
  onSelectPointGroupOverride,
  highlightedIrrepId = null,
  onHoverIrrep,
  onSelectMo,
}) => {
  const [showHeaderOptions, setShowHeaderOptions] = useState<boolean>(false);
  const [customChars, setCustomChars] = useState<number[]>(
    pointGroup.defaultReduciblePresets[0]?.characters ||
      pointGroup.classes.map((_, i) => (i === 0 ? 4 : 0))
  );
  const [activePresetLabel, setActivePresetLabel] = useState<string>(
    pointGroup.defaultReduciblePresets[0]?.label || 'Custom Γ_red'
  );

  // Sync default reducible characters when pointGroup changes
  useEffect(() => {
    const firstPreset = pointGroup.defaultReduciblePresets[0];
    if (firstPreset) {
      setCustomChars(firstPreset.characters);
      setActivePresetLabel(firstPreset.label);
    } else {
      setCustomChars(pointGroup.classes.map((_, i) => (i === 0 ? 4 : 0)));
      setActivePresetLabel('Custom Γ_red');
    }
  }, [pointGroup]);

  const decomposition = decomposeReducibleRepresentation(pointGroup, customChars);

  const directSumLatex =
    decomposition.results
      .filter((r) => r.coefficient > 0)
      .map(
        (r) =>
          `${r.coefficient === 1 ? '' : r.coefficient}${irrepToLatex(r.irrepSymbol)}`
      )
      .join(' \\oplus ') || '0';

  const symmetryElementLabelToLatex = (label: string): string => {
    return label
      .replace(/^C(\d+)\(([^)]+)\)/, 'C_{$1}($2)')
      .replace(/^C∞\(([^)]+)\)/, 'C_\\infty($1)')
      .replace(/^C2'\(([^)]+)\)/, "C_2'($1)")
      .replace(/^C2''\(([^)]+)\)/, "C_2''($1)")
      .replace(/^σv'\(([^)]+)\)/, "\\sigma_v'($1)")
      .replace(/^σv\(([^)]+)\)/, '\\sigma_v($1)')
      .replace(/^σh\(([^)]+)\)/, '\\sigma_h($1)')
      .replace(/^σd\(([^)]+)\)/, '\\sigma_d($1)')
      .replace(/^σ\(([^)]+)\)/, '\\sigma($1)');
  };

  return (
    <div className="bg-white border border-slate-200 rounded-lg flex flex-col overflow-hidden">
      {/* Header: Point Group Identity & Collapsible Symmetry Options */}
      <div className="px-4 py-3 bg-slate-50 border-b border-slate-200 flex flex-wrap items-center justify-between gap-3">
        <div>
          <div className="flex items-center gap-2">
            <h2 className="text-sm font-semibold text-slate-900">
              Symmetry Point Group &amp; Irreducible Representations
            </h2>
            <span className="text-slate-300" aria-hidden="true">·</span>
            <span className="font-mono text-xs font-bold text-blue-700">
              <Latex math={schoenfliesToLatex(pointGroup.schoenflies)} /> (Order{' '}
              <Latex
                math={`h = ${pointGroup.order >= 64 ? '\\infty' : pointGroup.order}`}
              />
              )
            </span>
          </div>
          <p className="text-xs text-slate-500 mt-0.5">{pointGroup.description}</p>
        </div>

        <button
          type="button"
          onClick={() => setShowHeaderOptions((v) => !v)}
          aria-expanded={showHeaderOptions}
          className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors flex items-center gap-1 whitespace-nowrap shrink-0 ${
            showHeaderOptions || highlightedSymmetryElementId
              ? 'bg-slate-900 text-white border-slate-900'
              : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-100'
          }`}
        >
          <SlidersHorizontal className="w-3 h-3" />
          <span>3D Symmetry Elements &amp; Group Selector</span>
          {showHeaderOptions ? (
            <ChevronUp className="w-3.5 h-3.5" />
          ) : (
            <ChevronDown className="w-3.5 h-3.5" />
          )}
        </button>
      </div>

      {/* Collapsible Interactive 3D Symmetry Elements & Group Selector Bar (Collapsed by Default) */}
      {showHeaderOptions && (
        <div className="px-4 py-2.5 bg-slate-100/70 border-b border-slate-200 flex flex-wrap items-center justify-between gap-3 text-xs">
          <div className="flex flex-wrap items-center gap-1.5">
            <span className="text-slate-600 font-medium mr-1">
              3D Symmetry Elements (click to highlight in 3D):
            </span>
            {pointGroup.elements3D.map((el) => {
              const isHi = highlightedSymmetryElementId === el.id;
              const parts = el.label.split(' ');
              const symToken = parts[0];
              const restLabel = parts.slice(1).join(' ');
              return (
                <button
                  key={el.id}
                  type="button"
                  onClick={() => onHighlightSymmetryElement(isHi ? null : el.id)}
                  className={`px-2.5 py-1 text-xs font-mono rounded border transition-colors whitespace-nowrap shrink-0 ${
                    isHi
                      ? 'bg-blue-600 text-white border-blue-600'
                      : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-50'
                  }`}
                >
                  <Latex math={symmetryElementLabelToLatex(symToken)} />
                  {restLabel ? ` ${restLabel}` : ''}
                </button>
              );
            })}
          </div>

          {onSelectPointGroupOverride && (
            <div className="flex items-center gap-2">
              <label
                htmlFor="pg-select"
                className="text-xs font-medium text-slate-600 whitespace-nowrap"
              >
                Inspect Group:
              </label>
              <select
                id="pg-select"
                value={pointGroup.id}
                onChange={(e) => onSelectPointGroupOverride(e.target.value)}
                className="px-2.5 py-1 text-xs font-mono font-semibold bg-white border border-slate-200 rounded text-slate-900 focus:outline-none focus:border-blue-600"
              >
                {Object.values(POINT_GROUPS).map((pg) => (
                  <option key={pg.id} value={pg.id}>
                    {pg.schoenflies} — {pg.name} (h={pg.order})
                  </option>
                ))}
              </select>
            </div>
          )}
        </div>
      )}

      {/* Character Table */}
      <div className="p-4 overflow-x-auto">
        <div className="flex items-center justify-between mb-2">
          <h3 className="text-xs font-semibold text-slate-800">
            Character Table for Point Group{' '}
            <Latex math={schoenfliesToLatex(pointGroup.schoenflies)} />
          </h3>
          <span className="text-[11px] text-slate-500">
            Highlighted row matches active orbital{' '}
            {selectedMo ? (
              <>
                <Latex math={moLabelToLatex(selectedMo.label)} /> (
                <Latex math={irrepToLatex(selectedMo.irrepSymbol)} />)
              </>
            ) : (
              '—'
            )}
          </span>
        </div>

        <table className="w-full text-left border-collapse border border-slate-200 text-xs font-mono tabular-nums">
          <thead>
            <tr className="bg-slate-100 text-slate-700 border-b border-slate-200">
              <th className="py-2 px-3 font-bold border-r border-slate-200">
                <Latex math={schoenfliesToLatex(pointGroup.schoenflies)} />
              </th>
              {pointGroup.classes.map((cls) => (
                <th
                  key={cls.id}
                  title={cls.description}
                  className="py-2 px-2.5 text-center font-semibold border-r border-slate-200"
                >
                  <Latex math={classSymbolToLatex(cls.symbol)} />
                </th>
              ))}
              <th className="py-2 px-3 font-sans font-semibold border-r border-slate-200">
                Linear / Rotations (<Latex math="x,y,z; R_\alpha" />)
              </th>
              <th className="py-2 px-3 font-sans font-semibold border-r border-slate-200">
                Quadratic / <Latex math="d" />-Orbitals
              </th>
              <th className="py-2 px-3 font-sans font-semibold">
                Spectroscopic Activity
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-slate-200">
            {pointGroup.irreps.map((irrep) => {
              const isSelectedMoIrrep =
                selectedMo &&
                (selectedMo.irrepId === irrep.id ||
                  selectedMo.irrepSymbol === irrep.symbol);
              const isHoveredIrrep =
                highlightedIrrepId &&
                (highlightedIrrepId === irrep.id ||
                  highlightedIrrepId === irrep.symbol);
              const matchingMos = molecule.molecularOrbitals.filter(
                (m) => m.irrepId === irrep.id || m.irrepSymbol === irrep.symbol
              );

              const linearStr =
                irrep.linearFunctions.length > 0
                  ? irrep.linearFunctions.join(', ')
                  : '—';
              const quadStr =
                irrep.orbitalBasis || irrep.quadraticFunctions.join(', ') || '—';

              return (
                <tr
                  key={irrep.id}
                  onMouseEnter={() => onHoverIrrep?.(irrep.id)}
                  onMouseLeave={() => onHoverIrrep?.(null)}
                  onClick={() => {
                    if (matchingMos.length > 0 && onSelectMo) {
                      onSelectMo(matchingMos[matchingMos.length - 1]);
                    }
                  }}
                  title={
                    matchingMos.length > 0
                      ? `Click to inspect ${matchingMos.map((m) => m.label).join(', ')} in 3D`
                      : `${irrep.symbol} symmetry species`
                  }
                  className={`${
                    matchingMos.length > 0 ? 'cursor-pointer' : ''
                  } transition-colors ${
                    isSelectedMoIrrep
                      ? 'bg-blue-50/90 text-slate-900 font-semibold'
                      : isHoveredIrrep
                        ? 'bg-teal-50/80 text-slate-900'
                        : 'bg-white text-slate-700 hover:bg-slate-50'
                  }`}
                >
                  <td className="py-2 px-3 font-bold border-r border-slate-200 text-blue-700 whitespace-nowrap">
                    <Latex math={irrepToLatex(irrep.symbol)} />
                    {isSelectedMoIrrep ? ' (Active MO)' : matchingMos.length > 0 ? ` (${matchingMos.length} MO)` : ''}
                  </td>
                  {irrep.characters.map((ch, idx) => (
                    <td
                      key={idx}
                      className={`py-2 px-2.5 text-center border-r border-slate-200 ${
                        ch < 0 ? 'text-rose-700' : 'text-slate-900'
                      }`}
                    >
                      {ch >= 0 ? `+${ch}` : `−${Math.abs(ch)}`}
                    </td>
                  ))}
                  <td className="py-2 px-3 border-r border-slate-200 text-slate-700">
                    {linearStr === '—' ? (
                      '—'
                    ) : (
                      <Latex math={basisFunctionToLatex(linearStr)} />
                    )}
                  </td>
                  <td className="py-2 px-3 border-r border-slate-200 text-slate-700">
                    {quadStr === '—' ? (
                      '—'
                    ) : (
                      <Latex math={basisFunctionToLatex(quadStr)} />
                    )}
                  </td>
                  <td className="py-2 px-3 font-sans text-[11px] whitespace-nowrap">
                    <span
                      className={
                        irrep.irActive
                          ? 'text-teal-700 font-semibold'
                          : 'text-slate-400'
                      }
                    >
                      {irrep.irActive ? 'IR Active' : 'IR Inactive'}
                    </span>
                    <span className="mx-1.5 text-slate-300" aria-hidden="true">
                      ·
                    </span>
                    <span
                      className={
                        irrep.ramanActive
                          ? 'text-blue-700 font-semibold'
                          : 'text-slate-400'
                      }
                    >
                      {irrep.ramanActive ? 'Raman Active' : 'Raman Inactive'}
                    </span>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Great Orthogonality Theorem: Reducible Representation (Γ_red) Decomposer */}
      <div className="p-4 bg-slate-50 border-t border-slate-200 space-y-3">
        <div className="flex flex-wrap items-center justify-between gap-2">
          <div>
            <h3 className="text-xs font-semibold text-slate-900">
              Reducible Representation (<Latex math="\Gamma_{\mathrm{red}}" />) Decomposition Calculator
            </h3>
            <p className="text-[11px] text-slate-600 mt-0.5">
              Great Orthogonality Theorem:{' '}
              <Latex math="a_i = \frac{1}{h} \sum_{R \in G} N_R \, \chi_{\mathrm{red}}(R) \, \chi_i^*(R)" />
            </p>
          </div>

          <div className="flex flex-wrap items-center gap-1.5">
            {pointGroup.defaultReduciblePresets.map((preset) => (
              <button
                key={preset.id}
                type="button"
                onClick={() => {
                  setCustomChars(preset.characters);
                  setActivePresetLabel(preset.label);
                }}
                className={`px-2.5 py-1 text-xs font-medium rounded border transition-colors whitespace-nowrap shrink-0 ${
                  activePresetLabel === preset.label
                    ? 'bg-slate-900 text-white border-slate-900'
                    : 'bg-white text-slate-700 border-slate-200 hover:bg-slate-100'
                }`}
              >
                <Latex math={reducibleLabelToLatex(preset.label.split(' ')[0])} />{' '}
                {preset.label.split(' ').slice(1).join(' ')}
              </button>
            ))}
          </div>
        </div>

        {/* Editable Character Inputs per Conjugacy Class */}
        <div className="grid grid-cols-2 sm:grid-cols-4 md:grid-cols-6 lg:grid-cols-8 gap-2">
          {pointGroup.classes.map((cls, idx) => (
            <div
              key={cls.id}
              className="bg-white border border-slate-200 rounded p-2 flex flex-col"
            >
              <label
                htmlFor={`char-input-${cls.id}`}
                className="text-[11px] font-mono font-semibold text-slate-600 mb-1"
              >
                <Latex math={`\\chi\\!\\left(${classSymbolToLatex(cls.symbol)}\\right)`} />{' '}
                <span className="text-slate-400">
                  [<Latex math={`N_R=${cls.count}`} />]
                </span>
              </label>
              <input
                id={`char-input-${cls.id}`}
                type="number"
                value={customChars[idx] ?? 0}
                onChange={(e) => {
                  const val = parseFloat(e.target.value) || 0;
                  const next = [...customChars];
                  next[idx] = val;
                  setCustomChars(next);
                  setActivePresetLabel('Custom Γ_red');
                }}
                className="w-full px-2 py-1 text-xs font-mono font-bold text-slate-900 bg-slate-50 border border-slate-200 rounded focus:outline-none focus:border-blue-600 tabular-nums"
              />
            </div>
          ))}
        </div>

        {/* Direct Sum Result & Step-by-Step Breakdown */}
        <div className="p-3 bg-white border border-slate-200 rounded-md flex flex-col sm:flex-row sm:items-center justify-between gap-3">
          <div>
            <span className="text-[11px] text-slate-500 block">
              Irreducible Representation Direct Sum ({activePresetLabel}):
            </span>
            <div
              className="text-base font-bold text-blue-700 mt-0.5"
              data-testid="irrep-decomposition-formula"
            >
              <Latex math={`\\Gamma_{\\mathrm{red}} = ${directSumLatex}`} />
            </div>
          </div>

          <div className="flex flex-wrap items-center gap-2 text-xs font-mono text-slate-700">
            {decomposition.results
              .filter((r) => r.coefficient > 0)
              .map((r) => (
                <div
                  key={r.irrepId}
                  className="bg-slate-50 px-2.5 py-1 rounded border border-slate-200"
                >
                  <Latex
                    math={`a\\!\\left(${irrepToLatex(r.irrepSymbol)}\\right) = \\frac{${r.sumContribution}}{${pointGroup.order}} = \\mathbf{${r.coefficient}}`}
                  />
                </div>
              ))}
          </div>
        </div>
      </div>
    </div>
  );
};
