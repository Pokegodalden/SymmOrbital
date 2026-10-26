import { AtomNode, BondEdge, ELEMENT_PROPERTIES, ElementSymbol } from './moleculesData';

export async function fetchWithRetry(
  url: string,
  options: RequestInit = {},
  maxRetries = 3,
  baseDelayMs = 600
): Promise<Response> {
  for (let attempt = 0; attempt <= maxRetries; attempt++) {
    const response = await fetch(url, options);
    if (response.ok) {
      return response;
    }
    if (
      (response.status === 429 || response.status >= 500) &&
      attempt < maxRetries
    ) {
      const retryAfter = Number(response.headers.get('Retry-After'));
      const delay =
        !Number.isNaN(retryAfter) && retryAfter > 0
          ? retryAfter * 1000
          : baseDelayMs * Math.pow(2, attempt);
      await new Promise((resolve) => setTimeout(resolve, delay));
      continue;
    }
    const errorBody = await response.text().catch(() => '');
    throw new Error(
      `HTTP ${response.status} (${response.statusText}) from ${url}: ${errorBody.slice(0, 200)}`
    );
  }
  throw new Error(`Exhausted retries for ${url}`);
}

export interface PubChemCompoundProfile {
  cid: number;
  title: string;
  iupacName: string;
  molecularFormula: string;
  molecularWeight: number;
  canonicalSmiles: string;
  inchikey: string;
  xlogp?: number;
  tpsa?: number;
  charge?: number;
  hBondDonorCount?: number;
  hBondAcceptorCount?: number;
  rotatableBondCount?: number;
  heavyAtomCount?: number;
  complexity?: number;
  conformerType: '3D-SDF' | '2D-SDF' | 'None';
  meanBondLengthAngstrom?: number;
  experimentalNotes: string[];
  image2dUrl: string;
  pubchemRecordUrl: string;
  atoms: AtomNode[];
  bonds: BondEdge[];
}

const SUPPORTED_ELEMENTS = new Set<string>(Object.keys(ELEMENT_PROPERTIES));

/**
 * Parses a PubChem V2000 SDF block into AtomNode[] and BondEdge[]
 * for our 2D Sketcher and 3D MO Engine.
 */
export function parseSdfToGraph(sdfText: string): {
  atoms: AtomNode[];
  bonds: BondEdge[];
  meanBondLengthAngstrom?: number;
} {
  const lines = sdfText.split(/\r?\n/);
  const countsLineIdx = lines.findIndex((l) => l.includes('V2000'));
  if (countsLineIdx === -1) {
    return { atoms: [], bonds: [] };
  }

  const countsLine = lines[countsLineIdx];
  const numAtoms = parseInt(countsLine.slice(0, 3).trim(), 10) || 0;
  const numBonds = parseInt(countsLine.slice(3, 6).trim(), 10) || 0;

  const rawAtoms: AtomNode[] = [];
  const elemCounts: Record<string, number> = {};

  for (let i = 1; i <= numAtoms; i++) {
    const line = lines[countsLineIdx + i];
    if (!line) break;
    const x = parseFloat(line.slice(0, 10).trim()) || 0;
    const y = parseFloat(line.slice(10, 20).trim()) || 0;
    const z = parseFloat(line.slice(20, 30).trim()) || 0;
    const rawSym = line.slice(31, 34).trim();
    const element: ElementSymbol = SUPPORTED_ELEMENTS.has(rawSym)
      ? (rawSym as ElementSymbol)
      : 'C';

    elemCounts[element] = (elemCounts[element] || 0) + 1;
    rawAtoms.push({
      id: `${element}${elemCounts[element]}`,
      element,
      x2d: Number(x.toFixed(2)),
      y2d: Number((-y).toFixed(2)),
      x3d: Number(x.toFixed(2)),
      y3d: Number(y.toFixed(2)),
      z3d: Number(z.toFixed(2)),
    });
  }

  const rawBonds: BondEdge[] = [];
  let bondLengthSum = 0;
  let bondLengthCount = 0;

  for (let i = 1; i <= numBonds; i++) {
    const line = lines[countsLineIdx + numAtoms + i];
    if (!line) break;
    const a1 = parseInt(line.slice(0, 3).trim(), 10) - 1;
    const a2 = parseInt(line.slice(3, 6).trim(), 10) - 1;
    const orderNum = parseInt(line.slice(6, 9).trim(), 10) || 1;
    const order: 1 | 2 | 3 = orderNum === 2 ? 2 : orderNum === 3 ? 3 : 1;

    if (rawAtoms[a1] && rawAtoms[a2]) {
      rawBonds.push({
        id: `sdf_b_${i}`,
        source: rawAtoms[a1].id,
        target: rawAtoms[a2].id,
        order,
      });
      const d = Math.hypot(
        rawAtoms[a1].x3d - rawAtoms[a2].x3d,
        rawAtoms[a1].y3d - rawAtoms[a2].y3d,
        rawAtoms[a1].z3d - rawAtoms[a2].z3d
      );
      if (d > 0.4) {
        bondLengthSum += d;
        bondLengthCount++;
      }
    }
  }

  const meanBondLengthAngstrom =
    bondLengthCount > 0 ? Number((bondLengthSum / bondLengthCount).toFixed(3)) : undefined;

  // Keep up to 18 atoms so the MO diagram & 2D canvas remain crisp and responsive
  if (rawAtoms.length > 18) {
    const nonH = rawAtoms.filter((a) => a.element !== 'H').slice(0, 14);
    const keepIds = new Set(nonH.map((a) => a.id));
    const filteredBonds = rawBonds.filter(
      (b) => keepIds.has(b.source) && keepIds.has(b.target)
    );
    return { atoms: nonH, bonds: filteredBonds, meanBondLengthAngstrom };
  }

  return { atoms: rawAtoms, bonds: rawBonds, meanBondLengthAngstrom };
}

/**
 * Extracts brief experimental property strings from a PubChem PUG-View JSON tree.
 */
function extractPugViewStrings(node: any, out: string[], maxCount = 4): void {
  if (!node || out.length >= maxCount) return;
  if (node.Value?.StringWithMarkup) {
    for (const item of node.Value.StringWithMarkup) {
      if (item?.String && typeof item.String === 'string' && item.String.length < 140) {
        out.push(item.String);
        if (out.length >= maxCount) return;
      }
    }
  }
  if (Array.isArray(node.Information)) {
    for (const info of node.Information) {
      extractPugViewStrings(info, out, maxCount);
      if (out.length >= maxCount) return;
    }
  }
  if (Array.isArray(node.Section)) {
    for (const sec of node.Section) {
      extractPugViewStrings(sec, out, maxCount);
      if (out.length >= maxCount) return;
    }
  }
}

/**
 * Resolves a chemical query (Name, CID, or SMILES) via live PubChem PUG REST & PUG-View APIs
 * and returns molecular properties, 2D/3D atomic coordinates, and experimental annotations.
 */
export async function lookupPubChemCompound(query: string): Promise<PubChemCompoundProfile> {
  const clean = query.trim();
  if (!clean) {
    throw new Error('Please enter a chemical name, formula, SMILES, or PubChem CID.');
  }

  const isNumericCid = /^\d+$/.test(clean);
  const namespace = isNumericCid ? 'cid' : 'name';
  const encoded = encodeURIComponent(clean);

  const propFields =
    'Title,IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,ConnectivitySMILES,InChIKey,XLogP,TPSA,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,Complexity';

  // 1. Fetch compound properties from PubChem PUG REST
  const propUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${namespace}/${encoded}/property/${propFields}/JSON`;
  let propRes: Response;
  try {
    propRes = await fetchWithRetry(propUrl);
  } catch {
    // Fallback to SMILES namespace if name lookup failed
    const smilesUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encoded}/property/${propFields}/JSON`;
    propRes = await fetchWithRetry(smilesUrl);
  }

  const propData = await propRes.json();
  const prop = propData?.PropertyTable?.Properties?.[0];
  if (!prop || !prop.CID) {
    throw new Error(`No PubChem record found for "${clean}".`);
  }

  const cid = Number(prop.CID);

  // 2. Fetch 3D Conformer SDF (fallback to 2D SDF if 3D conformer unavailable)
  let sdfText = '';
  let conformerType: '3D-SDF' | '2D-SDF' | 'None' = 'None';
  try {
    const sdf3dRes = await fetchWithRetry(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`,
      {},
      1
    );
    sdfText = await sdf3dRes.text();
    conformerType = '3D-SDF';
  } catch {
    try {
      const sdf2dRes = await fetchWithRetry(
        `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=2d`,
        {},
        1
      );
      sdfText = await sdf2dRes.text();
      conformerType = '2D-SDF';
    } catch {
      sdfText = '';
      conformerType = 'None';
    }
  }

  // 3. Optionally fetch PUG-View Experimental Properties summary
  const experimentalNotes: string[] = [];
  try {
    const pugViewRes = await fetchWithRetry(
      `https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/${cid}/JSON?heading=Experimental+Properties`,
      {},
      0
    );
    const pugViewJson = await pugViewRes.json();
    extractPugViewStrings(pugViewJson?.Record, experimentalNotes, 3);
  } catch {
    // Optional annotation endpoint; ignore if compound lacks Experimental Properties heading
  }

  const parsedSdf = sdfText
    ? parseSdfToGraph(sdfText)
    : { atoms: [], bonds: [], meanBondLengthAngstrom: undefined };

  return {
    cid,
    title: prop.Title || clean,
    iupacName: prop.IUPACName || prop.Title || clean,
    molecularFormula: prop.MolecularFormula || '',
    molecularWeight: Number(prop.MolecularWeight) || 0,
    canonicalSmiles: prop.CanonicalSMILES || prop.ConnectivitySMILES || '',
    inchikey: prop.InChIKey || '',
    xlogp: prop.XLogP !== undefined ? Number(prop.XLogP) : undefined,
    tpsa: prop.TPSA !== undefined ? Number(prop.TPSA) : undefined,
    charge: prop.Charge !== undefined ? Number(prop.Charge) : 0,
    hBondDonorCount: prop.HBondDonorCount !== undefined ? Number(prop.HBondDonorCount) : undefined,
    hBondAcceptorCount: prop.HBondAcceptorCount !== undefined ? Number(prop.HBondAcceptorCount) : undefined,
    rotatableBondCount: prop.RotatableBondCount !== undefined ? Number(prop.RotatableBondCount) : undefined,
    heavyAtomCount: prop.HeavyAtomCount !== undefined ? Number(prop.HeavyAtomCount) : undefined,
    complexity: prop.Complexity !== undefined ? Number(prop.Complexity) : undefined,
    conformerType,
    meanBondLengthAngstrom: parsedSdf.meanBondLengthAngstrom,
    experimentalNotes,
    image2dUrl: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG?image_size=large`,
    pubchemRecordUrl: `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`,
    atoms: parsedSdf.atoms,
    bonds: parsedSdf.bonds,
  };
}
