#!/usr/bin/env python3
"""
03c_metacyc_mapping.py — EC indices -> MetaCyc pathway indices (authenticated)
===============================================================================
Track 3 Host-Microbe Integration

Companion to 03_functional_profiles.py. Maps the genome-resolved EC-number
pool indices (data/interim/functional_prediction/kegg_out/ec_abundance.tsv)
onto MetaCyc pathways via the BioCyc web services.

AUTHENTICATION (required — BioCyc gates these queries behind a session):
  The script expects BIOCYC_USER / BIOCYC_PASSWORD environment variables
  (it also accepts plain USER / PASSWORD, as defined in biocyc.env).
  Load them WITHOUT printing them, e.g.:
      set -a; source ~/biocyc.env; set +a
      python src/asv/03c_metacyc_mapping.py ...
  The session cookie is kept in memory only. Credentials are never written
  to disk, logged, or included in any output file.

METHOD
------
  1. Establish a BioCyc session (POST /credentials/login/).
  2. For each profiled EC number, one BioVelo query returning MetaCyc
     pathways that contain a reaction with that EC number
     (query: [x:x<-meta^^pathways,y<-(reactions-of-pathway x),
              "EC" instringci y^ec-number], detail=low).
  3. Aggregate: MetaCyc pathway pool index = sum of EC pool indices of its
     profiled reactions. DESCRIPTIVE ONLY (n=1 pooled per group; no stats).
  4. All responses cached under the api cache dir for offline re-runs.

Rate limit: >=1 s between requests, per BioCyc's published guidance.

USAGE
-----
  python src/asv/03c_metacyc_mapping.py \
      --ec data/interim/functional_prediction/kegg_out/ec_abundance.tsv \
      --out-dir data/interim/functional_prediction \
      [--cache-dir data/interim/literature_tables/api_cache] \
      [--max-ecs N]   # for smoke tests
"""

import argparse
import hashlib
import http.cookiejar
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

USER_AGENT = "Mozilla/5.0 (research; track3_host_microbe)"
BIOCYC = "https://websvc.biocyc.org"
SLEEP_S = 1.05  # BioCyc asks for <=1 request/sec on average


class FileCache:
    def __init__(self, cache_dir: Path):
        self.dir = Path(cache_dir)
        self.dir.mkdir(parents=True, exist_ok=True)

    def _path(self, key: str) -> Path:
        return self.dir / (hashlib.sha256(key.encode()).hexdigest() + ".txt")

    def get(self, key: str) -> str | None:
        p = self._path(key)
        return p.read_text(encoding="utf-8", errors="replace") if p.exists() else None

    def put(self, key: str, text: str) -> None:
        self._path(key).write_text(text, encoding="utf-8")


class BioCycSession:
    """Authenticated session with an in-memory cookie jar and disk cache."""

    def __init__(self, user: str, password: str, cache: FileCache):
        self.cache = cache
        self.jar = http.cookiejar.CookieJar()
        opener = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(self.jar))
        opener.addheaders = [("User-Agent", USER_AGENT)]
        self.opener = opener
        # login (never cached)
        data = urllib.parse.urlencode({"email": user,
                                       "password": password}).encode()
        req = urllib.request.Request(f"{BIOCYC}/credentials/login/",
                                     data=data, method="POST")
        with self.opener.open(req, timeout=60) as resp:
            self.login_status = resp.status
        if not any(c.name in ("PTools-session", "userIdentifier")
                   for c in self.jar):
            raise RuntimeError(
                "BioCyc login did not establish a session cookie — check "
                "BIOCYC_USER / BIOCYC_PASSWORD.")

    def get(self, url: str, retries: int = 4, timeout: int = 120) -> str:
        cached = self.cache.get(url)
        if cached is not None:
            return cached
        time.sleep(SLEEP_S)
        last_err = None
        for attempt in range(retries):
            try:
                with self.opener.open(url, timeout=timeout) as resp:
                    text = resp.read().decode("utf-8", errors="replace")
                self.cache.put(url, text)
                return text
            except Exception as exc:  # noqa: BLE001
                last_err = exc
                time.sleep(3.0 * (attempt + 1))
        raise RuntimeError(f"BioCyc GET failed after {retries} tries: "
                           f"{url[:120]}... ({last_err})")


def query_pathways_for_ec(sess: BioCycSession, ec: str) -> list[tuple[str, str]]:
    """Return [(pathway_frame_id, pathway_common_name)] for one EC number."""
    q = (f'[x:x<-meta^^pathways,y<-(reactions-of-pathway x),'
         f'"{ec}" instringci y^ec-number]')
    url = (f"{BIOCYC}/xmlquery?query={urllib.parse.quote(q)}&detail=low")
    text = sess.get(url)
    out = []
    try:
        root = ET.fromstring(text)
        for el in root.findall(".//Pathway"):
            fid = el.get("frameid")
            if not fid:
                continue
            nm = ""
            cn = el.find("common-name")
            if cn is not None and cn.text:
                nm = cn.text.strip()
            out.append((fid, nm))
    except ET.ParseError:
        pass
    # dedupe frames, prefer entries carrying a name
    by_frame: dict[str, str] = {}
    for fid, nm in out:
        if fid not in by_frame or (nm and not by_frame[fid]):
            by_frame[fid] = nm
    return list(by_frame.items())


def fetch_all_pathway_names(sess: BioCycSession) -> dict[str, str]:
    """One request: all MetaCyc pathways at low detail -> frame -> name."""
    q = "[x:x<-meta^^pathways]"
    url = f"{BIOCYC}/xmlquery?query={urllib.parse.quote(q)}&detail=low"
    text = sess.get(url)
    names: dict[str, str] = {}
    try:
        root = ET.fromstring(text)
        for el in root.findall(".//Pathway"):
            fid = el.get("frameid")
            cn = el.find("common-name")
            if fid and cn is not None and cn.text:
                names[fid] = cn.text.strip()
    except ET.ParseError:
        pass
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ec", required=True,
                    help="EC abundance table from 03_functional_profiles.py")
    ap.add_argument("--out-dir", default="data/interim/functional_prediction")
    ap.add_argument("--cache-dir",
                    default="data/interim/literature_tables/api_cache")
    ap.add_argument("--max-ecs", type=int, default=0,
                    help="Limit ECs queried (smoke tests only)")
    args = ap.parse_args()

    user = (os.environ.get("BIOCYC_USER") or os.environ.get("USER") or "").strip()
    password = (os.environ.get("BIOCYC_PASSWORD")
                or os.environ.get("PASSWORD") or "").strip()
    if not user or not password:
        print("[META] ERROR: set BIOCYC_USER / BIOCYC_PASSWORD (or USER / "
              "PASSWORD, as in biocyc.env) and re-run, e.g. "
              "'set -a; source ~/biocyc.env; set +a'", file=sys.stderr)
        return 1

    base = Path(__file__).resolve().parents[2]
    ec_path = Path(args.ec)
    out_dir = Path(args.out_dir)
    if not ec_path.is_absolute():
        ec_path = base / ec_path
    if not out_dir.is_absolute():
        out_dir = base / out_dir
    metacyc_out = out_dir / "metacyc_out"
    metacyc_out.mkdir(parents=True, exist_ok=True)
    cache = FileCache(Path(args.cache_dir) if Path(args.cache_dir).is_absolute()
                      else base / args.cache_dir)

    ec_df = pd.read_csv(ec_path, sep="\t")
    ecs = [str(e).replace("ec:", "") for e in ec_df["ec"].tolist()]
    if args.max_ecs and args.max_ecs > 0:
        ecs = ecs[: args.max_ecs]
    print(f"[META] {len(ecs)} EC numbers to map (from {ec_path})")

    print("[META] Logging into BioCyc ...")
    sess = BioCycSession(user, password, cache)
    print(f"[META] Login HTTP {sess.login_status} — session established")

    ec_idx = {str(r["ec"]).replace("ec:", ""): (float(r["jumper_idx"]),
                                              float(r["laggard_idx"]))
              for _, r in ec_df.iterrows()}

    # EC -> pathways (one BioVelo query per EC, cached)
    pairs: dict[str, set[str]] = {}          # ec -> {pathway frames}
    pw_names: dict[str, str] = {}
    n_hits = 0
    for i, ec in enumerate(ecs, 1):
        try:
            hits = query_pathways_for_ec(sess, ec)
        except Exception as exc:  # noqa: BLE001
            print(f"[META] EC {ec} query failed: {exc}")
            continue
        if hits:
            pairs[ec] = {h[0] for h in hits}
            for h in hits:
                if h[1]:
                    pw_names.setdefault(h[0], h[1])
            n_hits += 1
        if i % 100 == 0 or i == len(ecs):
            print(f"[META]   {i}/{len(ecs)} ECs queried "
                  f"({n_hits} with MetaCyc hits)", flush=True)

    all_frames = {f for fs in pairs.values() for f in fs}
    missing_names = all_frames - set(pw_names)
    if missing_names:
        print(f"[META] {len(missing_names)} pathways lack names — fetching "
              f"the full MetaCyc pathway list (one cached request) ...")
        try:
            extra = fetch_all_pathway_names(sess)
            for f, nm in extra.items():
                if f in missing_names and nm:
                    pw_names[f] = nm
            still = sum(1 for f in missing_names if f not in pw_names)
            print(f"[META] names resolved for {len(missing_names) - still} "
                  f"pathways ({still} still unnamed)")
        except Exception as exc:  # noqa: BLE001
            print(f"[META] pathway-name fetch failed: {exc}")

    print(f"[META] {n_hits}/{len(ecs)} ECs map to MetaCyc pathways "
          f"({len(all_frames)} unique pathways)")

    # pathway indices: sum of EC indices per pathway
    pw_j: dict[str, float] = {}
    pw_l: dict[str, float] = {}
    pw_ec_count: dict[str, int] = {}
    for ec, frames in pairs.items():
        j, l = ec_idx.get(ec, (0.0, 0.0))
        for f in frames:
            pw_j[f] = pw_j.get(f, 0.0) + j
            pw_l[f] = pw_l.get(f, 0.0) + l
            pw_ec_count[f] = pw_ec_count.get(f, 0) + 1

    pc = 1e-6
    rows = []
    for f in sorted(pw_j):
        j, l = pw_j[f], pw_l[f]
        ratio = (j + pc) / (l + pc)
        import math
        rows.append({
            "metacyc_pathway": f,
            # class frames (e.g. Respiration) have no common-name:
            # fall back to the readable frame id
            "pathway_name": pw_names.get(f, "") or f,
            "n_ec_mapped": pw_ec_count[f],
            "jumper_idx": j,
            "laggard_idx": l,
            "log2_ratio_idx": math.log2(ratio),
        })
    out_df = pd.DataFrame(rows)
    if len(out_df):
        out_df = out_df.sort_values("log2_ratio_idx", key=abs, ascending=False)
    else:
        out_df = pd.DataFrame(columns=["metacyc_pathway", "pathway_name",
                                       "n_ec_mapped", "jumper_idx",
                                       "laggard_idx", "log2_ratio_idx"])
    out_path = metacyc_out / "metacyc_pathway_abundance.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[META] Wrote {len(out_df)} MetaCyc pathway indices -> {out_path}")

    (metacyc_out / "README.txt").write_text(
        "MetaCyc pathway indices — method & provenance\n"
        "=============================================\n"
        "Generated by src/asv/03c_metacyc_mapping.py (authenticated BioCyc\n"
        "web-services session; credentials read from BIOCYC_USER / \n"
        "BIOCYC_PASSWORD environment variables and never written to disk).\n\n"
        "Method: for each genome-resolved EC pool index (KEGG GENOME -> KO ->\n"
        "reaction -> EC chain, see kegg_out/ec_abundance.tsv), one BioVelo\n"
        "query retrieves the MetaCyc pathways containing a reaction with that\n"
        "EC number. Pathway pool index = sum of EC pool indices of its\n"
        "profiled reactions.\n\n"
        "DESCRIPTIVE ONLY: n=1 pooled sample per group. No p-values, FDR or\n"
        "statistical inference were computed from the two pool columns.\n\n"
        f"Run stats: {n_hits}/{len(ecs)} ECs mapped; {len(out_df)} pathways "
        f"profiled. Pathways without a common-name are MetaCyc class frames "
        f"(e.g. Respiration, Fermentation) and use their frame id as name.
"
        "All responses cached under data/interim/literature_tables/api_cache/.
")

    print("[META] DONE — no p-values were computed (n=1 design).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
