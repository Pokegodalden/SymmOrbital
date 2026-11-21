/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import React, { useEffect, useRef, useState } from 'react';
import { Check, ChevronDown, Layers } from 'lucide-react';

export type AppVersion = 'v3' | 'v2' | 'v1';

interface VersionOption {
  id: AppVersion;
  shortLabel: string;
  title: string;
  subtitle: string;
  tag: string;
}

const VERSION_OPTIONS: VersionOption[] = [
  {
    id: 'v3',
    shortLabel: 'v3.0',
    title: 'v3.0 · Computational Workbench',
    subtitle: 'Löwdin EHT, 27 groups, MEP/slice, Walsh & XYZ/Cube I/O',
    tag: 'Latest',
  },
  {
    id: 'v2',
    shortLabel: 'v2.0',
    title: 'v2.0 · Pro 3-Pane Workbench',
    subtitle: 'Synchronized 3-pane stage, Character Table & NIST audit',
    tag: 'Previous',
  },
  {
    id: 'v1',
    shortLabel: 'v1.0',
    title: 'v1.0 · Classic Stacked Layout',
    subtitle: 'Original vertical single-column workflow snapshot',
    tag: 'Legacy',
  },
];

interface VersionSwitcherDropdownProps {
  currentVersion: AppVersion;
  onSelectVersion: (version: AppVersion) => void;
}

export function VersionSwitcherDropdown({
  currentVersion,
  onSelectVersion,
}: VersionSwitcherDropdownProps) {
  const [open, setOpen] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);

  const activeOption =
    VERSION_OPTIONS.find((opt) => opt.id === currentVersion) || VERSION_OPTIONS[0];

  useEffect(() => {
    if (!open) return;

    const handlePointerDown = (e: MouseEvent) => {
      if (containerRef.current && !containerRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };

    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        setOpen(false);
      }
    };

    document.addEventListener('mousedown', handlePointerDown);
    document.addEventListener('keydown', handleKeyDown);
    return () => {
      document.removeEventListener('mousedown', handlePointerDown);
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [open]);

  return (
    <div ref={containerRef} className="relative shrink-0">
      <button
        type="button"
        onClick={() => setOpen((prev) => !prev)}
        aria-haspopup="menu"
        aria-expanded={open}
        aria-label={`Switch SymmOrbital version (current: ${activeOption.shortLabel})`}
        data-testid="version-switcher-button"
        className="px-2.5 py-1.5 text-xs font-medium text-slate-700 bg-slate-100 border border-slate-200 rounded-md hover:bg-slate-200 transition-colors flex items-center gap-1.5 whitespace-nowrap shrink-0"
      >
        <Layers className="w-3.5 h-3.5 text-blue-600 shrink-0" />
        <span className="font-mono font-semibold">{activeOption.shortLabel}</span>
        <ChevronDown
          className={`w-3.5 h-3.5 text-slate-500 transition-transform ${
            open ? 'rotate-180' : ''
          }`}
        />
      </button>

      {open && (
        <div
          role="menu"
          aria-label="SymmOrbital Versions"
          data-testid="version-switcher-menu"
          className="absolute right-0 mt-1.5 w-80 bg-white border border-slate-200 rounded-lg shadow-lg py-1.5 z-50"
        >
          <div className="px-3 py-1.5 border-b border-slate-100 flex items-center justify-between text-[11px] font-mono text-slate-500">
            <span>SYMMORBITAL VERSIONS</span>
            <span>Switch anytime</span>
          </div>

          <div className="py-1">
            {VERSION_OPTIONS.map((option) => {
              const isSelected = option.id === currentVersion;
              return (
                <button
                  key={option.id}
                  type="button"
                  role="menuitemradio"
                  aria-checked={isSelected}
                  data-testid={`version-option-${option.id}`}
                  onClick={() => {
                    onSelectVersion(option.id);
                    setOpen(false);
                  }}
                  className={`w-full text-left px-3 py-2 flex items-start gap-2.5 transition-colors ${
                    isSelected
                      ? 'bg-blue-50 text-slate-900'
                      : 'text-slate-700 hover:bg-slate-50'
                  }`}
                >
                  <div className="w-4 h-4 mt-0.5 flex items-center justify-center shrink-0">
                    {isSelected && <Check className="w-3.5 h-3.5 text-blue-600" />}
                  </div>
                  <div className="min-w-0 flex-1">
                    <div className="flex items-center justify-between gap-2">
                      <span className="text-xs font-semibold text-slate-900 truncate">
                        {option.title}
                      </span>
                      <span
                        className={`text-[10px] font-mono shrink-0 ${
                          option.id === 'v3'
                            ? 'text-blue-700 font-semibold'
                            : 'text-slate-500'
                        }`}
                      >
                        {option.tag}
                      </span>
                    </div>
                    <p className="text-[11px] text-slate-500 mt-0.5 leading-snug">
                      {option.subtitle}
                    </p>
                  </div>
                </button>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}
