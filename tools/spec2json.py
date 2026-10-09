"""Convert the extracted spec text files into per-suite JSON used by the app."""
import sys, re, json, os

src, out = sys.argv[1], sys.argv[2]
JUNK = {"DNV", "IEC", "+1", "+3", "Part 1: Data schema", "Open Geospatial Consortium"}


def lines(name):
    with open(os.path.join(src, name), encoding="utf8") as f:
        return [l.strip() for l in f.read().split("\n") if l.strip() and l.strip() not in JUNK]


suites = {n: {"n": n} for n in range(1, 18)}
common = {}

# ---- modules (analyses) -------------------------------------------------
cur = None
for l in lines("modules.offshore.txt"):
    m = re.match(r"^(\d{1,2})\. (.+)$", l)
    if m and 1 <= int(m.group(1)) <= 16:
        cur = suites[int(m.group(1))]
        cur["title"] = m.group(2)
        cur["analyses"] = []
        continue
    if cur is None:
        continue
    if l.startswith("Applicable to"):
        cur["scope"] = l
    else:
        cur["analyses"].append(l.rstrip("."))
suites[17]["title"] = "Engineering Economics and Techno-Economic Analysis"

# ---- equations / models -------------------------------------------------
cur, sect = None, None
uni = []
for l in lines("equations.models.txt"):
    m = re.match(r"^Suite (\d{1,2})\. (.+)$", l)
    if m:
        cur = suites[int(m.group(1))]
        cur["models"] = []
        continue
    if l.startswith("Universal computational framework"):
        cur = "uni"
        continue
    if cur is None:
        continue
    if cur == "uni":
        uni.append(l)
        continue
    m = re.match(r"^([A-Z][^.]{3,90})\. (.+)$", l)
    if m:
        cur["models"].append({"group": m.group(1), "text": m.group(2)})
    else:
        cur["models"].append({"group": "", "text": l})
common["framework"] = uni

# ---- inputs / outputs ---------------------------------------------------
cur, mode = None, None
tail = []
for l in lines("input.output.offshore.txt"):
    m = re.match(r"^Suite (\d{1,2})\. (.+)$", l)
    if m:
        cur = suites[int(m.group(1))]
        cur["inputs"], cur["outputs"] = [], []
        mode = "sub"
        continue
    if l.startswith("Integration of the 17 Suites"):
        cur = "tail"
    if cur is None:
        continue
    if cur == "tail":
        tail.append(l)
        continue
    if l.startswith("A. "):
        mode = "in"
        continue
    if l.startswith("B. "):
        mode = "out"
        continue
    if mode == "sub":
        cur["subtitle"] = l
    elif mode == "in":
        m = re.match(r"^([^:]{3,60}): (.+)$", l)
        if m:
            cur["inputs"].append({"cat": m.group(1), "text": m.group(2)})
        else:
            cur.setdefault("ioNote", l)
    elif mode == "out":
        cur["outputs"].append(l.rstrip("."))
common["integration"] = tail

# ---- initial / boundary conditions --------------------------------------
cur, mode = None, None
tail = []
for l in lines("ibcs.offshore.txt"):
    m = re.match(r"^Suite (\d{1,2})\. (.+)$", l)
    if m:
        cur = suites[int(m.group(1))]
        cur["ic"], cur["bc"] = [], []
        mode = None
        continue
    if l.startswith("Common Initial-Condition"):
        cur = "tail"
    if cur is None:
        continue
    if cur == "tail":
        tail.append(l)
        continue
    if l.startswith("A. "):
        mode = "ic"
        continue
    if l.startswith("B. "):
        mode = "bc"
        continue
    if mode:
        cur[mode].append(l)
common["conditions"] = tail

# ---- calibration / verification / validation ----------------------------
cur, mode = None, None
tail = []
for l in lines("cvv.offshore.txt"):
    m = re.match(r"^Suite (\d{1,2})\. (.+)$", l)
    if m:
        cur = suites[int(m.group(1))]
        cur["cvv"] = {"calibration": [], "verification": [], "validation": []}
        mode = None
        continue
    if l.startswith("Common Calibration, Verification"):
        cur = "tail"
    if cur is None:
        continue
    if cur == "tail":
        tail.append(l)
        continue
    if l in ("Calibration", "Verification", "Validation"):
        mode = l.lower()
        continue
    if mode:
        cur["cvv"][mode].append(l)
common["cvv"] = tail

# ---- geometry formats ---------------------------------------------------
cur = None
tail = []
for l in lines("geometry.formats.offshore.txt"):
    m = re.match(r"^(\d{1,2})\. (.+)$", l)
    if m and 1 <= int(m.group(1)) <= 16:
        cur = suites[int(m.group(1))]
        cur["geometry"] = []
        continue
    if l.startswith("Essential universal geometry-import"):
        cur = "tail"
    if cur is None:
        continue
    if cur == "tail":
        tail.append(l)
        continue
    m = re.match(r"^([^:]{3,80}): (.+)$", l)
    if m:
        cur["geometry"].append({"cat": m.group(1), "text": m.group(2)})
common["geometry"] = tail

suites[17]["analyses"] = [m["group"].rstrip(".") + " analysis" for m in suites[17]["models"] if m["group"]]
suites[17]["scope"] = "Operates across the other 16 suites: receives quantities, schedules, energy yield, reliability and design life and converts them into financial and economic evaluations."
os.makedirs(out, exist_ok=True)
index = []
for n, s in suites.items():
    s.setdefault("analyses", [])
    s.setdefault("geometry", [])
    with open(os.path.join(out, f"suite-{n:02d}.json"), "w", encoding="utf8") as f:
        json.dump(s, f, ensure_ascii=False, separators=(",", ":"))
    index.append({"n": n, "title": s["title"], "subtitle": s.get("subtitle", ""), "scope": s.get("scope", ""),
                  "analyses": s["analyses"]})
    print(n, s["title"], "| analyses", len(s["analyses"]), "models", len(s.get("models", [])), "in", len(s.get("inputs", [])),
          "out", len(s.get("outputs", [])), "ic", len(s.get("ic", [])), "bc", len(s.get("bc", [])),
          "cvv", [len(v) for v in s.get("cvv", {}).values()], "geo", len(s["geometry"]))
with open(os.path.join(out, "index.json"), "w", encoding="utf8") as f:
    json.dump(index, f, ensure_ascii=False, separators=(",", ":"))
with open(os.path.join(out, "common.json"), "w", encoding="utf8") as f:
    json.dump(common, f, ensure_ascii=False, separators=(",", ":"))
print({k: len(v) for k, v in common.items()})
