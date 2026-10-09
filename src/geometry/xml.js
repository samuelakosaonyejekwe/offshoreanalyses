// Small tolerant XML reader producing a plain object tree. DOMParser does not exist in Web Workers or Node,
// and the geospatial / civil formats read here (LandXML, GML, KML, GPX, 3MF, VTK XML, SVG) only need elements,
// attributes and text. Namespaces are not resolved: every element keeps its qualified `tag` and its
// prefix-less local `name`, which is what the importers match on.
//
// Not supported (by design): DTD validation, external entities, custom entity declarations (ignored).

const NAMED = { lt: '<', gt: '>', amp: '&', quot: '"', apos: "'", nbsp: ' ' };

/** Replace the predefined and numeric character references in `s`. Unknown entities are left as written. */
export function decodeEntities(s) {
  if (s.indexOf('&') < 0) return s;
  return s.replace(/&(#x[0-9a-fA-F]+|#\d+|[A-Za-z][A-Za-z0-9]*);/g, (m, e) => {
    if (e[0] === '#') {
      const code = e[1] === 'x' || e[1] === 'X' ? parseInt(e.slice(2), 16) : parseInt(e.slice(1), 10);
      return Number.isFinite(code) && code > 0 && code <= 0x10ffff ? String.fromCodePoint(code) : m;
    }
    return Object.prototype.hasOwnProperty.call(NAMED, e) ? NAMED[e] : m;
  });
}

function localOf(tag) {
  const c = tag.indexOf(':');
  return c >= 0 ? tag.slice(c + 1) : tag;
}

/**
 * Parse XML text.
 * @param {string} text
 * @param {object} [opts]
 * @param {boolean} [opts.keepWhitespace=false] keep whitespace-only text between elements
 * @returns {{ tag: string, name: string, attrs: Object<string,string>, children: object[], text: string }}
 *   the root element. `text` is the concatenation of the element's own character data (CDATA included,
 *   entities decoded, not trimmed). Attribute keys are the qualified names as written; use `attr()` to look
 *   one up by local name.
 * @throws {Error} when no element is found
 */
export function parseXml(text, opts = {}) {
  const keepWs = !!opts.keepWhitespace;
  const n = text.length;
  const doc = { tag: '#document', name: '#document', attrs: {}, children: [], text: '' };
  const stack = [doc];
  let cur = doc;
  let pos = 0;
  const addText = (raw, decode) => {
    if (cur === doc) return;
    if (!keepWs && !/\S/.test(raw)) return;
    cur.text += decode ? decodeEntities(raw) : raw;
  };
  while (pos < n) {
    const lt = text.indexOf('<', pos);
    if (lt < 0) { addText(text.slice(pos), true); break; }
    if (lt > pos) addText(text.slice(pos, lt), true);
    const c = text.charCodeAt(lt + 1);
    if (c === 33) { // '!'
      if (text.startsWith('<!--', lt)) { const e = text.indexOf('-->', lt + 4); pos = e < 0 ? n : e + 3; continue; }
      if (text.startsWith('<![CDATA[', lt)) { const e = text.indexOf(']]>', lt + 9); addText(text.slice(lt + 9, e < 0 ? n : e), false); pos = e < 0 ? n : e + 3; continue; }
      // DOCTYPE and other declarations, possibly with an internal subset in [ ... ]
      let depth = 0, i = lt + 2;
      for (; i < n; i++) { const ch = text.charCodeAt(i); if (ch === 91) depth++; else if (ch === 93) depth--; else if (ch === 62 && depth <= 0) break; }
      pos = i + 1;
      continue;
    }
    if (c === 63) { const e = text.indexOf('?>', lt + 2); pos = e < 0 ? n : e + 2; continue; } // processing instruction
    if (c === 47) { // closing tag
      const e = text.indexOf('>', lt + 2);
      const tag = text.slice(lt + 2, e < 0 ? n : e).trim();
      for (let k = stack.length - 1; k > 0; k--) {
        if (stack[k].tag === tag) { stack.length = k; cur = stack[k - 1]; break; }
      }
      pos = e < 0 ? n : e + 1;
      continue;
    }
    // opening tag: scan the name, then attributes, honouring quotes
    let i = lt + 1;
    while (i < n) { const ch = text.charCodeAt(i); if (ch === 32 || ch === 9 || ch === 10 || ch === 13 || ch === 62 || ch === 47) break; i++; }
    const tag = text.slice(lt + 1, i);
    if (!tag) { addText('<', false); pos = lt + 1; continue; }
    const el = { tag, name: localOf(tag), attrs: {}, children: [], text: '' };
    let selfClose = false;
    for (;;) {
      while (i < n) { const ch = text.charCodeAt(i); if (ch === 32 || ch === 9 || ch === 10 || ch === 13) i++; else break; }
      if (i >= n) break;
      const ch = text.charCodeAt(i);
      if (ch === 62) { i++; break; }
      if (ch === 47) { selfClose = true; i++; continue; }
      let j = i;
      while (j < n) { const cj = text.charCodeAt(j); if (cj === 61 || cj === 32 || cj === 9 || cj === 10 || cj === 13 || cj === 62 || cj === 47) break; j++; }
      const key = text.slice(i, j);
      while (j < n && (text.charCodeAt(j) === 32 || text.charCodeAt(j) === 9 || text.charCodeAt(j) === 10 || text.charCodeAt(j) === 13)) j++;
      if (text.charCodeAt(j) === 61) {
        j++;
        while (j < n && (text.charCodeAt(j) === 32 || text.charCodeAt(j) === 9 || text.charCodeAt(j) === 10 || text.charCodeAt(j) === 13)) j++;
        const q = text.charCodeAt(j);
        if (q === 34 || q === 39) {
          const e = text.indexOf(q === 34 ? '"' : "'", j + 1);
          const end = e < 0 ? n : e;
          el.attrs[key] = decodeEntities(text.slice(j + 1, end));
          i = end + 1;
        } else {
          let e = j;
          while (e < n) { const ce = text.charCodeAt(e); if (ce === 32 || ce === 9 || ce === 10 || ce === 13 || ce === 62) break; e++; }
          el.attrs[key] = decodeEntities(text.slice(j, e));
          i = e;
        }
      } else {
        if (key) el.attrs[key] = '';
        i = j === i ? j + 1 : j;
      }
    }
    cur.children.push(el);
    if (!selfClose) { stack.push(el); cur = el; }
    pos = i;
  }
  const roots = doc.children;
  if (!roots.length) throw new Error("Not an XML document: expected an element start tag such as '<Root ...>' but found none.");
  return roots.length === 1 ? roots[0] : { tag: '#document', name: '#document', attrs: {}, children: roots, text: '' };
}

/** Attribute by local name (prefix ignored), or `dflt`. */
export function attr(node, name, dflt) {
  const a = node.attrs;
  if (name in a) return a[name];
  for (const k in a) { const c = k.indexOf(':'); if (c >= 0 && k.slice(c + 1) === name) return a[k]; }
  return dflt;
}

/** Numeric attribute by local name, or `dflt` when absent or not a number. */
export function numAttr(node, name, dflt) {
  const v = attr(node, name);
  if (v === undefined || v === '') return dflt;
  const x = parseFloat(v);
  return Number.isFinite(x) ? x : dflt;
}

/** Direct children with the given local name (all element children when `name` is omitted). */
export function kids(node, name) {
  if (!name) return node.children;
  const out = [];
  for (const c of node.children) if (c.name === name) out.push(c);
  return out;
}

/** First direct child with the given local name, or null. */
export function kid(node, name) {
  for (const c of node.children) if (c.name === name) return c;
  return null;
}

/** All descendants (document order, excluding `node`) whose local name is `name` or in the array `name`. */
export function desc(node, name) {
  const want = Array.isArray(name) ? new Set(name) : null;
  const out = [];
  const stack = [];
  for (let i = node.children.length - 1; i >= 0; i--) stack.push(node.children[i]);
  while (stack.length) {
    const c = stack.pop();
    if (want ? want.has(c.name) : c.name === name) out.push(c);
    for (let i = c.children.length - 1; i >= 0; i--) stack.push(c.children[i]);
  }
  return out;
}

/** First descendant with the given local name, or null. */
export function firstDesc(node, name) {
  const stack = [];
  for (let i = node.children.length - 1; i >= 0; i--) stack.push(node.children[i]);
  while (stack.length) {
    const c = stack.pop();
    if (c.name === name) return c;
    for (let i = c.children.length - 1; i >= 0; i--) stack.push(c.children[i]);
  }
  return null;
}

/** Follow a path of local names through direct children: path(node, 'a', 'b') → node/a/b or null. */
export function path(node, ...names) {
  let c = node;
  for (const nm of names) { c = c && kid(c, nm); if (!c) return null; }
  return c;
}

/** Own text plus the text of all descendants, trimmed. */
export function deepText(node) {
  if (!node) return '';
  if (!node.children.length) return node.text.trim();
  let s = node.text;
  const stack = [...node.children].reverse();
  while (stack.length) { const c = stack.pop(); s += ' ' + c.text; for (let i = c.children.length - 1; i >= 0; i--) stack.push(c.children[i]); }
  return s.trim();
}

/** Trimmed text of the first direct child called `name`, or `dflt`. */
export function kidText(node, name, dflt) {
  const c = kid(node, name);
  return c ? c.text.trim() : dflt;
}

/**
 * Parse a whitespace- and/or comma-separated list of numbers without building an intermediate token array
 * (coordinate lists in GML / LandXML / VTK XML can hold millions of values).
 * @returns {Float64Array}
 */
export function numberList(s) {
  const n = s.length;
  let out = new Float64Array(Math.max(8, n >> 3));
  let count = 0, i = 0;
  while (i < n) {
    let c = s.charCodeAt(i);
    while (i < n && (c === 32 || c === 44 || c === 10 || c === 13 || c === 9 || c === 59)) c = s.charCodeAt(++i);
    if (i >= n) break;
    const st = i;
    while (i < n) { c = s.charCodeAt(i); if (c === 32 || c === 44 || c === 10 || c === 13 || c === 9 || c === 59) break; i++; }
    if (count === out.length) { const b = new Float64Array(out.length * 2); b.set(out); out = b; }
    out[count++] = +s.slice(st, i);
  }
  return out.slice(0, count);
}
