// DWG reading — worker side.
//
// This file is the whole of what the application contributes to reading a DWG drawing inside the worker: a
// message loop. On a "read" request that carries the user's consent it imports the DWG library from the address
// it is given (the library's own public distribution, see dwg.js), hands it the bytes of the drawing, and posts
// back what the library returns as plain data: header variables, layer names, block definitions and entities.
// No code of the library is part of this repository or of the application's modules; it is fetched by the
// browser, runs in this worker only, and exchanges nothing with the application but these messages.
//
// Request   { id, cmd: 'read', consent: true, url, bytes: ArrayBuffer }
// Response  { id, ok: true, drawing: { header, layers, blocks, entities } }  |  { id, ok: false, error }

self.onmessage = async (ev) => {
  const m = ev.data;
  if (!m || m.cmd !== 'read') return;
  try {
    if (m.consent !== true) throw new Error('the request does not carry the consent of the user');
    const lib = await import(m.url);
    const reader = await lib.LibreDwg.create();
    const ptr = reader.dwg_read_data(m.bytes, lib.Dwg_File_Type.DWG);
    if (ptr === undefined || ptr === null) throw new Error('the library could not open the data as a DWG drawing');
    let db;
    try { db = reader.convert(ptr); } finally { reader.dwg_free(ptr); }
    const header = {};
    for (const k in db.header || {}) { const v = db.header[k]; if (typeof v === 'number' || typeof v === 'string') header[k] = v; }
    const tables = db.tables || {};
    const plain = (list) => JSON.parse(JSON.stringify(list || [], (k, v) => (v instanceof Uint8Array || typeof v === 'bigint' ? undefined : v)));
    self.postMessage({
      id: m.id, ok: true,
      drawing: {
        header,
        layers: ((tables.LAYER && tables.LAYER.entries) || []).map((l) => l.name),
        blocks: ((tables.BLOCK_RECORD && tables.BLOCK_RECORD.entries) || []).map((b) => ({ name: b.name, basePoint: b.basePoint, entities: plain(b.entities) })),
        entities: plain(db.entities),
      },
    });
  } catch (e) {
    self.postMessage({ id: m.id, ok: false, error: String((e && e.message) || e) });
  }
};
