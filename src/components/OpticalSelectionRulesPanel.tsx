/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useMemo, useState } from 'react';
import { AtomNode, MolecularOrbital } from '../chemistry/moleculesData';
import { PointGroupData } from '../chemistry/symmetryData';
import {
  computeOpticalSelectionRules,
  evaluateElectronicTransition,
  wavelengthToVisibleColor,
} from '../chemistry/symmetryV3';
import {
  irrepToLatex,
  moLabelToLatex,
  schoenfliesToLatex,
} from '../chemistry/latexUtils';
import { Latex } from './Latex';
import { Sparkles } from 'lucide-react';

interface OpticalSelectionRulesPanelProps {
  pointGroup: PointGroupData;
  molecularOrbitals: MolecularOrbital[];
  atoms: AtomNode[];
  selectedMo: MolecularOrbital;
  onSelectMo: (mo: MolecularOrbital) => void;
}

export const OpticalSelectionRulesPanel: React.FC<OpticalSelectionRulesPanelProps> = ({
  pointGroup,
  molecularOrbitals,
  atoms,
  selectedMo,
  onSelectMo,
}) => {
  const occupied = useMemo(
    () => molecularOrbitals.filter((m) => m.occupancy > 0).sort((a, b) => b.energy - a.energy),
    [molecularOrbitals]
  );
  const virtuals = useMemo(
    () => molecularOrbitals.filter((m) => m.occupancy < 2).sort((a, b) => a.energy - b.energy),
    [molecularOrbitals]
  );

  const defaultFromId =
    occupied.find((m) => m.isHOMO)?.id || occupied[0]?.id || molecularOrbitals[0]?.id || '';
  const defaultToId =
    virtuals.find((m) => m.isLUMO)?.id ||
    virtuals.find((m) => m.id !== defaultFromId)?.id ||
    molecularOrbitals[molecularOrbitals.length - 1]?.id ||
    '';

  const [fromMoId, setFromMoId] = useState<string>(defaultFromId);
  const [toMoId, setToMoId] = useState<string>(defaultToId);

  const activeFromMo = useMemo(() => {
    return (
      molecularOrbitals.find((m) => m.id === fromMoId) ||
      molecularOrbitals.find((m) => m.id === defaultFromId) ||
      molecularOrbitals[0]
    );
  }, [molecularOrbitals, fromMoId, defaultFromId]);

  const activeToMo = useMemo(() => {
    const found = molecularOrbitals.find((m) => m.id === toMoId);
    if (found && found.id !== activeFromMo?.id) return found;
    return (
      molecularOrbitals.find((m) => m.id === defaultToId) ||
      molecularOrbitals[molecularOrbitals.length - 1]
    );
  }, [molecularOrbitals, toMoId, defaultToId, activeFromMo?.id]);

  const customRule = useMemo(() => {
    if (!activeFromMo || !activeToMo) return null;
    return evaluateElectronicTransition(pointGroup, activeFromMo, activeToMo, atoms);
  }, [pointGroup, activeFromMo, activeToMo, atoms]);

  const summary = useMemo(
    () => computeOpticalSelectionRules(pointGroup, molecularOrbitals, atoms),
    [pointGroup, molecularOrbitals, atoms]
  );

  return (
    <div
      data-testid="optical-selection-rules-panel"
      className="bg-white border border-slate-200 rounded-lg p-4 space-y-4"
    >
      {/* Header */}
      <div className="flex flex-wrap items-center justify-between gap-2 pb-3 border-b border-slate-200">
        <div>
          <h2 className="text-sm font-semibold text-slate-900 flex items-center gap-2">
            <Sparkles className="w-4 h-4 text-blue-600 shrink-0" />
            <span>
              Optical Selection Rules · UV-Vis Electronic Polarization, Laporte Parity &amp; IR / Raman Activity (
              <Latex math={schoenfliesToLatex(pointGroup.schoenflies)} />)
            </span>
          </h2>
          <p className="text-xs text-slate-500 mt-0.5">
            Exact direct-product transition dipole integral check{' '}
            <Latex math="\Gamma_i \otimes \Gamma_{\hat{\mu}} \otimes \Gamma_f \supset \Gamma_1" />
          </p>
        </div>

        <div className="text-xs font-mono text-slate-600">
          {summary.vibrational.isCentrosymmetric ? (
            <span className="text-blue-700 font-semibold">
              Centrosymmetric · Mutual Exclusion Active
            </span>
          ) : (
            <span>Non-centrosymmetric · IR/Raman Coincidence Allowed</span>
          )}
        </div>
      </div>

      {/* Interactive Transition Evaluator + Vibrational Summary */}
      <div className="grid grid-cols-1 lg:grid-cols-12 gap-4">
        {/* UV-Vis Pair Evaluator */}
        <div className="lg:col-span-7 border border-slate-200 rounded-lg p-3.5 space-y-3">
          <div className="flex flex-wrap items-center justify-between gap-2">
            <span className="text-xs font-semibold text-slate-900">
              Interactive Electronic Transition Calculator
            </span>
            <span className="text-[11px] font-mono text-slate-500">
              Active 3D MO: <Latex math={moLabelToLatex(selectedMo.label)} />
            </span>
          </div>

          <div className="grid grid-cols-1 sm:grid-cols-2 gap-2.5">
            <label className="text-xs space-y-1 block">
              <span className="text-[11px] font-medium text-slate-600 block">
                Initial Orbital (<Latex math="\psi_i" />)
              </span>
              <select
                aria-label="Initial molecular orbital"
                data-testid="optical-from-mo"
                value={activeFromMo?.id || ''}
                onChange={(e) => setFromMoId(e.target.value)}
                className="w-full px-2.5 py-1.5 text-xs font-mono bg-slate-50 border border-slate-300 rounded text-slate-900"
              >
                {occupied.map((mo) => (
                  <option key={mo.id} value={mo.id}>
                    {mo.label} ({mo.irrepSymbol}) · {mo.energy.toFixed(2)} eV
                    {mo.isHOMO ? ' [HOMO]' : ''}
                  </option>
                ))}
              </select>
            </label>

            <label className="text-xs space-y-1 block">
              <span className="text-[11px] font-medium text-slate-600 block">
                Final Orbital (<Latex math="\psi_f" />)
              </span>
              <select
                aria-label="Final molecular orbital"
                data-testid="optical-to-mo"
                value={activeToMo?.id || ''}
                onChange={(e) => setToMoId(e.target.value)}
                className="w-full px-2.5 py-1.5 text-xs font-mono bg-slate-50 border border-slate-300 rounded text-slate-900"
              >
                {virtuals.map((mo) => (
                  <option key={mo.id} value={mo.id}>
                    {mo.label} ({mo.irrepSymbol}) · {mo.energy.toFixed(2)} eV
                    {mo.isLUMO ? ' [LUMO]' : ''}
                  </option>
                ))}
              </select>
            </label>
          </div>

          {customRule && (
            <div className="bg-slate-50 border border-slate-200 rounded-md p-3 space-y-2.5 text-xs">
              <div className="flex flex-wrap items-center justify-between gap-2">
                <div className="font-mono text-slate-900 font-semibold">
                  <Latex math={moLabelToLatex(customRule.fromLabel)} /> (
                  <Latex math={irrepToLatex(customRule.fromIrrepSymbol)} />) →{' '}
                  <Latex math={moLabelToLatex(customRule.toLabel)} /> (
                  <Latex math={irrepToLatex(customRule.toIrrepSymbol)} />)
                </div>

                <div
                  data-testid="optical-verdict"
                  className={`font-mono font-semibold text-xs ${
                    customRule.isElectricDipoleAllowed
                      ? 'text-teal-700'
                      : 'text-amber-700'
                  }`}
                >
                  {customRule.isElectricDipoleAllowed
                    ? `Electric-Dipole Allowed (${customRule.allowedPolarizations.join(', ')}-polarized)`
                    : 'Electric-Dipole Forbidden (Vibronically Weak)'}
                </div>
              </div>

              <div className="flex flex-wrap items-center gap-2 text-[11px] font-mono text-slate-700 tabular-nums">
                <span>
                  Direct Product:{' '}
                  <strong className="text-slate-900">{customRule.directProductFormula}</strong>
                </span>
                <span className="text-slate-300" aria-hidden="true">·</span>
                <span>
                  <Latex math={`\\Delta E = ${customRule.deltaEv.toFixed(2)}\\;\\mathrm{eV}`} />
                </span>
                {customRule.wavelengthNm > 0 && (() => {
                  const customVis = wavelengthToVisibleColor(customRule.wavelengthNm);
                  return (
                    <>
                      <span className="text-slate-300" aria-hidden="true">·</span>
                      <span className="inline-flex items-center gap-1.5">
                        {customVis && (
                          <span
                            className="inline-block w-2.5 h-2.5 rounded-xs border border-black/15 shrink-0"
                            style={{ backgroundColor: customVis.rgb }}
                            title={`Visible ${customVis.bandName}`}
                          />
                        )}
                        <span>
                          <Latex math={`\\lambda = ${customRule.wavelengthNm.toFixed(1)}\\;\\mathrm{nm}`} /> (
                          {Math.round(customRule.wavenumberCm1).toLocaleString()} cm⁻¹
                          {customVis ? ` · ${customVis.bandName}` : ''})
                        </span>
                      </span>
                    </>
                  );
                })()}
                <span className="text-slate-300" aria-hidden="true">·</span>
                <span>
                  Laporte:{' '}
                  <strong
                    className={
                      customRule.laporteStatus === 'allowed'
                        ? 'text-teal-700'
                        : customRule.laporteStatus === 'forbidden'
                          ? 'text-amber-700'
                          : 'text-slate-600'
                    }
                  >
                    {customRule.laporteStatus === 'allowed'
                      ? 'Allowed (g ↔ u)'
                      : customRule.laporteStatus === 'forbidden'
                        ? 'Forbidden (same parity)'
                        : 'N/A (no i)'}
                  </strong>
                </span>
              </div>

              <p className="text-xs text-slate-600 leading-relaxed">{customRule.rationale}</p>

              <div className="flex items-center gap-2 pt-1">
                {activeFromMo && (
                  <button
                    type="button"
                    onClick={() => onSelectMo(activeFromMo)}
                    className="px-2.5 py-1 text-[11px] font-medium bg-white border border-slate-300 rounded hover:bg-slate-100 text-slate-700 transition-colors"
                  >
                    Inspect <Latex math={moLabelToLatex(activeFromMo.label)} /> in 3D
                  </button>
                )}
                {activeToMo && (
                  <button
                    type="button"
                    onClick={() => onSelectMo(activeToMo)}
                    className="px-2.5 py-1 text-[11px] font-medium bg-white border border-slate-300 rounded hover:bg-slate-100 text-slate-700 transition-colors"
                  >
                    Inspect <Latex math={moLabelToLatex(activeToMo.label)} /> in 3D
                  </button>
                )}
              </div>
            </div>
          )}
        </div>

        {/* IR & Raman Vibrational Selection Rules */}
        <div className="lg:col-span-5 border border-slate-200 rounded-lg p-3.5 space-y-2.5 flex flex-col justify-between">
          <div className="space-y-2.5">
            <h3 className="text-xs font-semibold text-slate-900">
              IR &amp; Raman Vibrational Symmetry Rules
            </h3>
            <p className="text-[11px] text-slate-500">
              IR-active modes transform as linear coordinates (<Latex math="x, y, z" />); Raman-active modes transform as quadratic polarizability tensors (<Latex math="x^2, xy, \dots" />).
            </p>

            <div className="space-y-2 text-xs font-mono">
              <div className="flex items-start justify-between gap-2 py-1.5 border-b border-slate-100">
                <span className="text-slate-500 shrink-0">IR Active:</span>
                <span className="text-right text-blue-700 font-semibold">
                  {summary.vibrational.irActiveIrreps.length > 0
                    ? summary.vibrational.irActiveIrreps.join(', ')
                    : 'None'}
                </span>
              </div>
              <div className="flex items-start justify-between gap-2 py-1.5 border-b border-slate-100">
                <span className="text-slate-500 shrink-0">Raman Active:</span>
                <span className="text-right text-teal-700 font-semibold">
                  {summary.vibrational.ramanActiveIrreps.length > 0
                    ? summary.vibrational.ramanActiveIrreps.join(', ')
                    : 'None'}
                </span>
              </div>
              <div className="flex items-start justify-between gap-2 py-1.5 border-b border-slate-100">
                <span className="text-slate-500 shrink-0">Both (Coincident):</span>
                <span className="text-right text-slate-800">
                  {summary.vibrational.coincidentIrreps.length > 0
                    ? summary.vibrational.coincidentIrreps.join(', ')
                    : 'None (Mutual Exclusion)'}
                </span>
              </div>
              <div className="flex items-start justify-between gap-2 py-1.5">
                <span className="text-slate-500 shrink-0">Silent (Inactive):</span>
                <span className="text-right text-slate-500">
                  {summary.vibrational.silentIrreps.length > 0
                    ? summary.vibrational.silentIrreps.join(', ')
                    : 'None'}
                </span>
              </div>
            </div>
          </div>

          <div className="text-[11px] text-slate-500 pt-2 border-t border-slate-100">
            {summary.vibrational.mutualExclusionApplies
              ? 'Because this point group contains inversion (i), no fundamental vibration can be simultaneously IR and Raman active.'
              : 'Lacking an inversion center (i), polar modes can appear in both IR absorption and Raman scattering spectra.'}
          </div>
        </div>
      </div>

      {/* Frontier Electronic Transitions Table */}
      {summary.frontierTransitions.length > 0 && (
        <div className="border border-slate-200 rounded-lg overflow-hidden">
          <div className="px-3 py-2 bg-slate-50 border-b border-slate-200 flex flex-wrap items-center justify-between gap-2">
            <span className="text-xs font-semibold text-slate-800">
              Lowest-Energy Frontier Electronic Transitions (Click row to evaluate)
            </span>
            <span className="text-[10.5px] font-mono text-slate-500">
              Visible band (380–750 nm) rows tinted by photon λ
            </span>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-xs font-mono tabular-nums">
              <thead className="bg-slate-50 text-slate-500 border-b border-slate-200">
                <tr>
                  <th className="text-left px-3 py-1.5 font-medium">Transition</th>
                  <th className="text-left px-2 py-1.5 font-medium">Direct Product</th>
                  <th className="text-right px-2 py-1.5 font-medium">ΔE (eV)</th>
                  <th className="text-right px-2 py-1.5 font-medium">λ (nm)</th>
                  <th className="text-left px-2 py-1.5 font-medium">Polarization</th>
                  <th className="text-left px-3 py-1.5 font-medium">Selection Verdict</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-slate-100">
                {summary.frontierTransitions.map((tr) => {
                  const isCurrent =
                    tr.fromMoId === activeFromMo?.id && tr.toMoId === activeToMo?.id;
                  const visColor = wavelengthToVisibleColor(tr.wavelengthNm);
                  const rowStyle = visColor
                    ? {
                        backgroundColor: `rgba(${visColor.r}, ${visColor.g}, ${visColor.b}, ${
                          isCurrent ? 0.34 : 0.22
                        })`,
                        boxShadow: isCurrent
                          ? `inset 4px 0 0 ${visColor.rgb}, inset 0 0 0 1px rgba(${visColor.r}, ${visColor.g}, ${visColor.b}, 0.72)`
                          : `inset 4px 0 0 ${visColor.rgb}`,
                      }
                    : undefined;
                  return (
                    <tr
                      key={`${tr.fromMoId}-${tr.toMoId}`}
                      data-visible-band={visColor ? visColor.bandName.toLowerCase() : undefined}
                      onClick={() => {
                        setFromMoId(tr.fromMoId);
                        setToMoId(tr.toMoId);
                      }}
                      style={rowStyle}
                      className={`cursor-pointer transition-colors ${
                        visColor
                          ? ''
                          : isCurrent
                            ? 'bg-blue-50'
                            : 'hover:bg-slate-50'
                      }`}
                    >
                      <td className="px-3 py-1.5 font-semibold text-slate-900">
                        <Latex math={moLabelToLatex(tr.fromLabel)} /> →{' '}
                        <Latex math={moLabelToLatex(tr.toLabel)} />
                      </td>
                      <td className="px-2 py-1.5 text-slate-600">{tr.directProductFormula}</td>
                      <td className="px-2 py-1.5 text-right text-slate-800">
                        {tr.deltaEv.toFixed(2)}
                      </td>
                      <td className="px-2 py-1.5 text-right text-slate-700">
                        {tr.wavelengthNm > 0 ? (
                          <span
                            className="inline-flex items-center justify-end gap-1.5"
                            title={
                              visColor
                                ? `Visible (${visColor.bandName}) · ${tr.wavelengthNm.toFixed(1)} nm`
                                : `${tr.wavelengthNm < 380 ? 'UV' : 'IR'} · ${tr.wavelengthNm.toFixed(1)} nm`
                            }
                          >
                            {visColor && (
                              <span
                                className="inline-block w-2.5 h-2.5 rounded-xs border border-black/15 shrink-0"
                                style={{ backgroundColor: visColor.rgb }}
                                aria-hidden="true"
                              />
                            )}
                            <span>{tr.wavelengthNm.toFixed(1)}</span>
                          </span>
                        ) : (
                          '—'
                        )}
                      </td>
                      <td className="px-2 py-1.5 text-blue-700 font-semibold">
                        {tr.allowedPolarizations.length > 0
                          ? tr.allowedPolarizations.join(', ')
                          : '—'}
                      </td>
                      <td
                        className={`px-3 py-1.5 font-semibold ${
                          tr.isElectricDipoleAllowed ? 'text-teal-700' : 'text-amber-700'
                        }`}
                      >
                        {tr.isElectricDipoleAllowed ? 'Allowed' : 'Forbidden'}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
};
