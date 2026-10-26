import { LCAOContribution } from './moleculesData';

/**
 * Converts a plain molecular formula (e.g. "H2O", "C6H6", "SF6", "CH2O")
 * into a LaTeX string: "\mathrm{H_2O}", "\mathrm{C_6H_6}".
 */
export function formulaToLatex(formula: string): string {
  if (!formula || formula === 'Empty') return '\\text{Empty}';
  const formatted = formula.replace(/([A-Za-z])(\d+)/g, '$1_{$2}');
  return `\\mathrm{${formatted}}`;
}

/**
 * Converts a Schoenflies point group symbol (e.g. "C2v", "D6h", "D∞h", "Td", "Oh", "Cs")
 * into standard LaTeX notation: "C_{2v}", "D_{6h}", "D_{\infty h}", "T_d", "O_h", "C_s".
 */
export function schoenfliesToLatex(pg: string): string {
  if (!pg) return '';
  const map: Record<string, string> = {
    C2v: 'C_{2v}',
    C3v: 'C_{3v}',
    D2h: 'D_{2h}',
    D3h: 'D_{3h}',
    D4h: 'D_{4h}',
    D6h: 'D_{6h}',
    Td: 'T_d',
    Oh: 'O_h',
    Cs: 'C_s',
    Cinfv: 'C_{\\infty v}',
    'C∞v': 'C_{\\infty v}',
    Dinfh: 'D_{\\infty h}',
    'D∞h': 'D_{\\infty h}',
  };
  if (map[pg]) return map[pg];
  // Generic fallback for Cnv / Dnh
  return pg.replace(/^([CDTSO])([0-9a-z∞]+)$/i, '$1_{$2}').replace('∞', '\\infty ');
}

/**
 * Converts a Mulliken Irreducible Representation symbol (e.g. "A1", "B2g", "E1u", "A2''", "Σg+", "Πu")
 * into standard LaTeX math notation: "A_1", "B_{2g}", "E_{1u}", "A_2''", "\Sigma_g^+", "\Pi_u".
 */
export function irrepToLatex(irrep: string): string {
  if (!irrep) return '';
  const trimmed = irrep.trim();

  // Direct sum e.g. "A1+B1+B2" or "A1 + B2"
  if (trimmed.includes('+')) {
    // Careful with Σ+ or Σg+ or Σu+
    if (trimmed === 'Σ+') return '\\Sigma^+';
    if (trimmed === 'Σg+') return '\\Sigma_g^+';
    if (trimmed === 'Σu+') return '\\Sigma_u^+';
    if (trimmed === '4Σ+ + 2Π') return '4\\Sigma^+ \\oplus 2\\Pi';
    if (trimmed === 'Σg+ + Σu+') return '\\Sigma_g^+ \\oplus \\Sigma_u^+';
    if (trimmed === 'Σu++Πu') return '\\Sigma_u^+ \\oplus \\Pi_u';
    if (trimmed === 'Σg++Σu++Πu+Πg') return '\\Sigma_g^+ \\oplus \\Sigma_u^+ \\oplus \\Pi_u \\oplus \\Pi_g';
    if (trimmed === 'Σ++Π') return '\\Sigma^+ \\oplus \\Pi';

    const parts = trimmed.split(/\s*\+\s*/);
    if (parts.length > 1) {
      return parts.map((p) => singleIrrepTermToLatex(p)).join(' \\oplus ');
    }
  }

  return singleIrrepTermToLatex(trimmed);
}

function singleIrrepTermToLatex(term: string): string {
  const m = term.trim().match(/^(\d*)([A-Za-zΣΠΔ]+)([0-9gu'+\-*]*)$/);
  if (!m) return term;
  const [, coeff, base, suffix] = m;

  let baseTex = base;
  if (base === 'Σ') baseTex = '\\Sigma';
  else if (base === 'Π') baseTex = '\\Pi';
  else if (base === 'Δ') baseTex = '\\Delta';

  if (!suffix) return `${coeff}${baseTex}`;

  // Handle primes and star superscripts vs subscripts
  let primes = '';
  let star = '';
  let plusMinus = '';
  let sub = suffix;

  if (sub.includes("''")) {
    primes = "''";
    sub = sub.replace("''", '');
  } else if (sub.includes("'")) {
    primes = "'";
    sub = sub.replace("'", '');
  }
  if (sub.includes('*')) {
    star = '*';
    sub = sub.replace(/\*/g, '');
  }
  if (sub.includes('+')) {
    plusMinus = '+';
    sub = sub.replace(/\+/g, '');
  } else if (sub.includes('-')) {
    plusMinus = '-';
    sub = sub.replace(/-/g, '');
  }

  const subTex = sub ? `_{${sub}}` : '';
  const supContent = `${primes}${plusMinus}${star}`;
  const supTex = supContent ? `^{${supContent}}` : '';

  return `${coeff}${baseTex}${subTex}${supTex}`;
}

/**
 * Converts a conjugacy class header symbol (e.g. "E", "2C3(z)", "3σv", "σv'(yz)", "2C∞(φ)", "3C2(=C4²)")
 * into standard LaTeX math notation.
 */
export function classSymbolToLatex(sym: string): string {
  if (!sym) return '';
  const exactMap: Record<string, string> = {
    E: 'E',
    i: 'i',
    'C2(z)': 'C_2(z)',
    'C2(y)': 'C_2(y)',
    'C2(x)': 'C_2(x)',
    'σv(xz)': '\\sigma_v(xz)',
    "σv'(yz)": "\\sigma_v'(yz)",
    '2C3(z)': '2C_3(z)',
    '3σv': '3\\sigma_v',
    'σ(xy)': '\\sigma(xy)',
    'σ(xz)': '\\sigma(xz)',
    'σ(yz)': '\\sigma(yz)',
    "3C2'": "3C_2'",
    "3C2''": "3C_2''",
    'σh(xy)': '\\sigma_h(xy)',
    '2S3': '2S_3',
    '2C4(z)': '2C_4(z)',
    "2C2'": "2C_2'",
    "2C2''": "2C_2''",
    '2S4': '2S_4',
    'σh': '\\sigma_h',
    '2σv': '2\\sigma_v',
    '2σd': '2\\sigma_d',
    '2C6(z)': '2C_6(z)',
    '2S6': '2S_6',
    '3σd': '3\\sigma_d',
    '8C3': '8C_3',
    '3C2': '3C_2',
    '6S4': '6S_4',
    '6σd': '6\\sigma_d',
    "6C2'": "6C_2'",
    '6C4': '6C_4',
    '3C2(=C4²)': '3C_2(=C_4^2)',
    '8S6': '8S_6',
    '3σh': '3\\sigma_h',
    '2C∞(φ)': '2C_\\infty(\\varphi)',
    '∞σv': '\\infty\\sigma_v',
    '2S∞(φ)': '2S_\\infty(\\varphi)',
    "∞C2'": "\\infty C_2'",
  };
  return exactMap[sym] || sym.replace(/σ/g, '\\sigma ').replace(/∞/g, '\\infty ');
}

/**
 * Converts a Molecular Orbital label (e.g. "1b1", "4a1*", "1e1g (π2)", "1e2u* (π4*)", "2b2 (n_O HOMO)")
 * into standard LaTeX notation: "1b_1", "4a_1^*", "1e_{1g}\\;(\\pi_2)", etc.
 */
export function moLabelToLatex(label: string): string {
  if (!label) return '';

  // Split optional parenthetical annotation e.g. "1e1g (π2)" -> main="1e1g", paren="π2"
  const m = label.trim().match(/^([^\s(]+)(?:\s*\((.+)\))?$/);
  if (!m) return label;

  const mainPart = formatSingleMoTokenLatex(m[1]);
  if (!m[2]) return mainPart;

  const parenRaw = m[2]
    .replace(/HOMO/g, '\\mathrm{HOMO}')
    .replace(/LUMO/g, '\\mathrm{LUMO}')
    .replace(/σ\*_([A-Za-z0-9]+)/g, '\\sigma^*_{\\mathrm{$1}}')
    .replace(/σ_([A-Za-z0-9]+)/g, '\\sigma_{\\mathrm{$1}}')
    .replace(/π\*_([A-Za-z0-9]+)/g, '\\pi^*_{\\mathrm{$1}}')
    .replace(/π_([A-Za-z0-9]+)/g, '\\pi_{\\mathrm{$1}}')
    .replace(/n_([A-Za-z0-9]+)/g, 'n_{\\mathrm{$1}}')
    .replace(/π(\d+)\*/g, '\\pi_{$1}^*')
    .replace(/π(\d+)/g, '\\pi_{$1}')
    .replace(/π\*/g, '\\pi^*')
    .replace(/σ\*/g, '\\sigma^*')
    .replace(/π/g, '\\pi')
    .replace(/σ/g, '\\sigma');

  return `${mainPart}\\;(${parenRaw})`;
}

function formatSingleMoTokenLatex(token: string): string {
  // e.g. "1b1", "4a1*", "1e1g", "1e2u*", "2a2''*", "4σg+", "3σu+", "1πu", "2πu*", "1e(a)"
  const subParen = token.match(/^([^()]+)\(([abxyz,]+)\)$/);
  if (subParen) {
    return `${formatSingleMoTokenLatex(subParen[1])}({${subParen[2]}})`;
  }

  const m = token.match(/^(\d+)(σ|π|δ|[a-zA-Z]+)([0-9gu'+\-*]*)$/);
  if (!m) return token;
  const [, num, sym, rest] = m;

  let symTex = sym;
  if (sym === 'σ') symTex = '\\sigma';
  else if (sym === 'π') symTex = '\\pi';
  else if (sym === 'δ') symTex = '\\delta';

  let primes = '';
  let star = '';
  let plus = '';
  let sub = rest;

  if (sub.includes("''")) {
    primes = "''";
    sub = sub.replace("''", '');
  } else if (sub.includes("'")) {
    primes = "'";
    sub = sub.replace("'", '');
  }
  if (sub.includes('*')) {
    star = '*';
    sub = sub.replace(/\*/g, '');
  }
  if (sub.includes('+')) {
    plus = '+';
    sub = sub.replace(/\+/g, '');
  }

  const subTex = sub ? `_{${sub}}` : '';
  const supContent = `${primes}${plus}${star}`;
  const supTex = supContent ? `^{${supContent}}` : '';

  return `${num}${symTex}${subTex}${supTex}`;
}

/**
 * Converts an LCAO wavefunction expansion into full LaTeX equation format:
 * \psi(1b_1) = 0.84\,\phi(\mathrm{O_1\;2s}) + 0.38\,\phi(\mathrm{H_1\;1s})
 */
export function lcaoExpansionToLatex(moLabel: string, lcao: LCAOContribution[]): string {
  const lhs = `\\psi\\!\\left(${moLabelToLatex(moLabel)}\\right) = `;
  if (!lcao || lcao.length === 0) return `${lhs} 0`;

  const terms = lcao.map((term, i) => {
    const sign = term.coeff >= 0 ? (i === 0 ? '' : ' + ') : i === 0 ? '-' : ' - ';
    const absCoeff = Math.abs(term.coeff).toFixed(2);
    const basisTex = aoLabelToLatex(term.label);
    return `${sign}${absCoeff}\\,\\phi\\!\\left(${basisTex}\\right)`;
  });

  return lhs + terms.join('');
}

/**
 * Converts an Atomic Orbital basis label (e.g. "O1 2pz", "C1 2pz(π)", "S 3dz²", "H1 1s")
 * into LaTeX notation.
 */
export function aoLabelToLatex(aoLabel: string): string {
  if (!aoLabel) return '';
  return aoLabel
    .replace(/^([A-Z][a-z]?)(\d+)\s+/, '\\mathrm{$1_{$2}}\\;')
    .replace(/^([A-Z][a-z]?)\s+/, '\\mathrm{$1}\\;')
    .replace(/(\d)p_?z/g, '$1p_z')
    .replace(/(\d)p_?x/g, '$1p_x')
    .replace(/(\d)p_?y/g, '$1p_y')
    .replace(/(\d)p_?σ/g, '$1p_\\sigma')
    .replace(/(\d)d_?z²|(\d)dz2/g, '$1$2d_{z^2}')
    .replace(/(\d)dx²-y²|(\d)dx2y2/g, '$1$2d_{x^2-y^2}')
    .replace(/\(in-plane\)/g, '\\;(\\text{in-plane})')
    .replace(/\(out-of-plane\)/g, '\\;(\\text{out-of-plane})')
    .replace(/\(π\)/g, '(\\pi)');
}

/**
 * Converts character table linear/rotational and quadratic/d-orbital basis strings into LaTeX.
 */
export function basisFunctionToLatex(basisStr: string): string {
  if (!basisStr || basisStr === '—') return '\\text{—}';
  return basisStr
    .replace(/Rx/g, 'R_x')
    .replace(/Ry/g, 'R_y')
    .replace(/Rz/g, 'R_z')
    .replace(/x²-y²/g, 'x^2-y^2')
    .replace(/x²\+y²\+z²/g, 'x^2+y^2+z^2')
    .replace(/2z²-x²-y²/g, '2z^2-x^2-y^2')
    .replace(/x²\+y²/g, 'x^2+y^2')
    .replace(/x²-3y²/g, 'x^2-3y^2')
    .replace(/x²/g, 'x^2')
    .replace(/y²/g, 'y^2')
    .replace(/z²/g, 'z^2')
    .replace(/dz²/g, 'd_{z^2}')
    .replace(/dx²-y²/g, 'd_{x^2-y^2}')
    .replace(/dxy/g, 'd_{xy}')
    .replace(/dxz/g, 'd_{xz}')
    .replace(/dyz/g, 'd_{yz}')
    .replace(/fxyz/g, 'f_{xyz}')
    .replace(/fx\(/g, 'f_x(')
    .replace(/fz\(/g, 'f_z(')
    .replace(/\bpx\b/g, 'p_x')
    .replace(/\bpy\b/g, 'p_y')
    .replace(/\bpz\b/g, 'p_z')
    .replace(/σg/g, '\\sigma_g')
    .replace(/σu/g, '\\sigma_u')
    .replace(/πg/g, '\\pi_g')
    .replace(/πu/g, '\\pi_u')
    .replace(/π\*/g, '\\pi^*')
    .replace(/π/g, '\\pi')
    .replace(/σ/g, '\\sigma')
    .replace(/δ/g, '\\delta')
    .replace(/—/g, '\\text{ — }');
}

/**
 * Converts a reducible representation preset label (e.g. "Γ_σ (2 O–H σ bonds)") into LaTeX.
 */
export function reducibleLabelToLatex(label: string): string {
  if (!label) return '';
  return label
    .replace(/Γ_val/g, '\\Gamma_{\\mathrm{val}}')
    .replace(/Γ_SALC/g, '\\Gamma_{\\mathrm{SALC}}')
    .replace(/Γ_vib/g, '\\Gamma_{\\mathrm{vib}}')
    .replace(/Γ_3N/g, '\\Gamma_{3N}')
    .replace(/Γ_σ/g, '\\Gamma_\\sigma')
    .replace(/Γ_π/g, '\\Gamma_\\pi')
    .replace(/Γ_CH/g, '\\Gamma_{\\mathrm{CH}}')
    .replace(/Γ_CC/g, '\\Gamma_{\\mathrm{CC}}')
    .replace(/Γ_red/g, '\\Gamma_{\\mathrm{red}}');
}
