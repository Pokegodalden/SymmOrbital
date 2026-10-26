import React, { useMemo } from 'react';
import katex from 'katex';

interface LatexProps {
  math: string;
  displayMode?: boolean;
  className?: string;
  title?: string;
}

/**
 * Renders chemical and mathematical expressions using KaTeX.
 */
export const Latex: React.FC<LatexProps> = ({
  math,
  displayMode = false,
  className = '',
  title,
}) => {
  const html = useMemo(() => {
    try {
      return katex.renderToString(math, {
        displayMode,
        throwOnError: false,
        strict: false,
        trust: true,
      });
    } catch {
      return math;
    }
  }, [math, displayMode]);

  return (
    <span
      className={`inline-block align-baseline ${className}`}
      title={title}
      dangerouslySetInnerHTML={{ __html: html }}
    />
  );
};

/**
 * Renders a molecular orbital or fragment label inside an SVG <text> element
 * with true mathematical subscript/superscript <tspan> shifts so SVG diagrams
 * match LaTeX typography without relying on external HTML overlays.
 */
export const SvgMathLabel: React.FC<{
  text: string;
  x: number;
  y: number;
  textAnchor?: 'start' | 'middle' | 'end';
  fill?: string;
  className?: string;
}> = ({ text, x, y, textAnchor = 'middle', fill = '#0F172A', className = '' }) => {
  // Parse tokens like "1b1", "4a1*", "1e1g (π2)", "1e2u* (π4*)", "2a1' (σ)", "3 × C 2pz (α)"
  const renderTokens = (raw: string) => {
    // Split into segments where digits/subscripts follow orbital letters
    const regex = /(\d+[a-zA-Zσπδ]+)([0-9gu'+\-*]+)|([πσ])(\d+\*?)|([A-Z][a-z]?)(\d+)|(2p|3d)(z²|x²-y²|[xyzσ])/g;
    const nodes: React.ReactNode[] = [];
    let lastIndex = 0;
    let match: RegExpExecArray | null;
    let key = 0;

    while ((match = regex.exec(raw)) !== null) {
      if (match.index > lastIndex) {
        nodes.push(raw.slice(lastIndex, match.index));
      }
      if (match[1] && match[2]) {
        // e.g. "1e" + "1g" or "4a" + "1*"
        const base = match[1];
        const suffix = match[2];
        const star = suffix.includes('*') ? '*' : '';
        const primes = suffix.includes("''") ? "''" : suffix.includes("'") ? "'" : '';
        const plus = suffix.includes('+') ? '+' : '';
        const sub = suffix.replace(/[*'+]/g, '');
        nodes.push(base);
        if (sub) {
          nodes.push(
            <tspan key={`sub_${key++}`} dy="2.2" fontSize="78%">
              {sub}
            </tspan>,
            <tspan key={`sub_reset_${key++}`} dy="-2.2">
              {'\u200B'}
            </tspan>
          );
        }
        const sup = `${primes}${plus}${star}`;
        if (sup) {
          nodes.push(
            <tspan key={`sup_${key++}`} dy="-3.2" fontSize="80%">
              {sup}
            </tspan>,
            <tspan key={`sup_reset_${key++}`} dy="3.2">
              {'\u200B'}
            </tspan>
          );
        }
      } else if (match[3] && match[4]) {
        // e.g. "π2" or "π4*"
        const sym = match[3];
        const subRaw = match[4];
        const hasStar = subRaw.endsWith('*');
        const digits = subRaw.replace('*', '');
        nodes.push(sym);
        nodes.push(
          <tspan key={`pisub_${key++}`} dy="2.2" fontSize="78%">
            {digits}
          </tspan>,
          <tspan key={`pisub_r_${key++}`} dy="-2.2">
            {'\u200B'}
          </tspan>
        );
        if (hasStar) {
          nodes.push(
            <tspan key={`pistar_${key++}`} dy="-3.2" fontSize="80%">
              *
            </tspan>,
            <tspan key={`pistar_r_${key++}`} dy="3.2">
              {'\u200B'}
            </tspan>
          );
        }
      } else if (match[5] && match[6]) {
        // e.g. "H1" or "C2"
        nodes.push(match[5]);
        nodes.push(
          <tspan key={`elsub_${key++}`} dy="2.2" fontSize="78%">
            {match[6]}
          </tspan>,
          <tspan key={`elsub_r_${key++}`} dy="-2.2">
            {'\u200B'}
          </tspan>
        );
      } else if (match[7] && match[8]) {
        // e.g. "2pz"
        nodes.push(match[7]);
        nodes.push(
          <tspan key={`aosub_${key++}`} dy="2.2" fontSize="78%">
            {match[8]}
          </tspan>,
          <tspan key={`aosub_r_${key++}`} dy="-2.2">
            {'\u200B'}
          </tspan>
        );
      }
      lastIndex = regex.lastIndex;
    }

    if (lastIndex < raw.length) {
      nodes.push(raw.slice(lastIndex));
    }
    return nodes;
  };

  return (
    <text x={x} y={y} textAnchor={textAnchor} fill={fill} className={className}>
      {renderTokens(text)}
    </text>
  );
};
