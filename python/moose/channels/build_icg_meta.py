#!/usr/bin/env python3
"""
build_icg_meta.py
=================
One-shot script that:

  1. Reads channel_db.csv.
  2. For each unique ModelDB modeldb_id, calls the ICGenealogy search API
     to discover the ICG internal channel IDs (one per channel suffix).
  3. Matches every (modeldb_id, suffix) row in channel_db.csv to its ICG ID.
  4. Fetches the ICG channel-detail record for each matched channel and
     collects the `ref` (paper metadata + citation count).
  5. Writes two output files:
       - channel_db.csv        updated in-place with a new `icg_id` column
       - icg_channel_meta.csv  one row per (icg_id): paper + citation data

The script is idempotent — rows that already have an icg_id are skipped
unless --force is given.  Re-run freely to refresh citation counts.

Usage
-----
    python build_icg_meta.py [--channel-db PATH] [--meta-csv PATH]
                             [--rate FLOAT] [--force] [--dry-run]

Defaults (relative to script location):
  --channel-db  data/channel_db.csv
  --meta-csv    data/icg_channel_meta.csv
  --rate        1.0   (requests per second)
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# ── API constants ─────────────────────────────────────────────────────────────

ICG_BASE   = 'https://icg.neurotheory.ox.ac.uk'
SEARCH_URL = ICG_BASE + '/api/app/search/?q={modeldb_id}'
DETAIL_URL = ICG_BASE + '/api/app/chs/{icg_id}/'

# ── default paths ─────────────────────────────────────────────────────────────

SCRIPT_DIR       = Path(__file__).parent
DEFAULT_CHAN_DB   = SCRIPT_DIR / 'data' / 'channel_db.csv'
DEFAULT_META_CSV  = SCRIPT_DIR / 'data' / 'icg_channel_meta.csv'

# ── output schema ─────────────────────────────────────────────────────────────

META_FIELDS = [
    'icg_id', 'modeldb_id', 'suffix', 'ion_class',
    'icg_ref_id', 'authors', 'title', 'year', 'pubmedid', 'citation_count',
]

# ── HTTP helpers ──────────────────────────────────────────────────────────────

def _fetch(url: str):
    req = Request(url, headers={'Accept': 'application/json'})
    try:
        with urlopen(req, timeout=20) as resp:
            return json.loads(resp.read().decode())
    except HTTPError as e:
        if e.code == 404:
            return None
        print(f'  [WARN] HTTP {e.code} for {url}', file=sys.stderr)
        return None
    except (URLError, OSError) as e:
        print(f'  [WARN] {e} for {url}', file=sys.stderr)
        return None


# ── ICG API helpers ───────────────────────────────────────────────────────────

def icg_search(modeldb_id: int, delay: float) -> list:
    """Return list of {icg_id, icg_suffix} for a ModelDB modeldb_id."""
    data = _fetch(SEARCH_URL.format(modeldb_id=modeldb_id))
    time.sleep(delay)
    if not isinstance(data, dict):
        return []

    prefix = str(modeldb_id) + '-'
    results = []
    for s in data.get('suggestions', []):
        val = s.get('value', '')
        if not val.startswith(prefix):
            continue
        # value: "87535-naxn [id: 2493, class: Na]"
        rest       = val[len(prefix):]           # "naxn [id: 2493, class: Na]"
        icg_suffix = rest.split(' ')[0]          # "naxn"
        try:
            icg_id = json.loads(s['data'])[0]    # [2493, 2]
        except Exception:
            continue
        results.append({'icg_id': icg_id, 'icg_suffix': icg_suffix})
    return results


def icg_detail(icg_id: int, delay: float) -> dict:
    """Return the `ref` dict from a channel detail, or {}."""
    data = _fetch(DETAIL_URL.format(icg_id=icg_id))
    time.sleep(delay)
    if not isinstance(data, dict):
        return {}
    return data.get('ref') or {}


# ── suffix matching ───────────────────────────────────────────────────────────

def _match(our: str, candidates: list) -> dict | None:
    """
    Find the best-matching ICG candidate for `our` suffix.
    Priority: exact > ours is prefix of theirs > theirs is prefix of ours
              > one contains the other.
    Returns the candidate dict or None.
    """
    lo = our.lower()
    for c in candidates:
        if c['icg_suffix'].lower() == lo:
            return c
    for c in candidates:
        if c['icg_suffix'].lower().startswith(lo):
            return c
    for c in candidates:
        if lo.startswith(c['icg_suffix'].lower()):
            return c
    for c in candidates:
        icg_lo = c['icg_suffix'].lower()
        if lo in icg_lo or icg_lo in lo:
            return c
    return None


# ── CSV helpers ───────────────────────────────────────────────────────────────

def load_channel_db(path: Path):
    """Return (fieldnames, rows) from channel_db.csv."""
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    return fieldnames, rows


def write_channel_db(path: Path, fieldnames: list, rows: list) -> None:
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def load_meta_csv(path: Path) -> dict:
    """Load existing icg_channel_meta.csv → {icg_id: row_dict}."""
    result = {}
    if not path.exists():
        return result
    with open(path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                result[int(row['icg_id'])] = row
            except (ValueError, KeyError):
                pass
    return result


def write_meta_csv(path: Path, meta: dict) -> None:
    rows = [meta[k] for k in sorted(meta)]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=META_FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--channel-db', default=str(DEFAULT_CHAN_DB))
    ap.add_argument('--meta-csv',   default=str(DEFAULT_META_CSV))
    ap.add_argument('--rate',  type=float, default=1.0,
                    help='API requests per second (default 1.0)')
    ap.add_argument('--force', action='store_true',
                    help='Re-fetch even rows that already have icg_id')
    ap.add_argument('--dry-run', action='store_true',
                    help='Fetch data but do not write any files')
    args = ap.parse_args()

    delay         = 1.0 / max(args.rate, 0.1)
    chan_db_path  = Path(args.channel_db)
    meta_csv_path = Path(args.meta_csv)

    fieldnames, rows = load_channel_db(chan_db_path)
    meta             = load_meta_csv(meta_csv_path)

    # Ensure icg_id column exists in fieldnames (insert after modeldb_id)
    if 'icg_id' not in fieldnames:
        idx = fieldnames.index('modeldb_id') + 1 if 'modeldb_id' in fieldnames else 0
        fieldnames = fieldnames[:idx] + ['icg_id'] + fieldnames[idx:]
        for r in rows:
            r.setdefault('icg_id', '')

    # Group rows by modeldb_id
    from collections import defaultdict
    by_model: dict[int, list] = defaultdict(list)
    for r in rows:
        try:
            by_model[int(r['modeldb_id'])].append(r)
        except (ValueError, TypeError):
            pass

    modeldb_ids = sorted(by_model)
    print(f'channel_db:  {len(rows)} rows, {len(modeldb_ids)} unique model IDs')
    print(f'meta cache:  {len(meta)} existing entries')
    print(f'Rate:        {args.rate:.1f} req/s (delay={delay:.2f}s)')
    print()

    matched_total = 0
    unmatched_total = 0
    new_detail_total = 0

    for idx, mid in enumerate(modeldb_ids, 1):
        model_rows = by_model[mid]

        # Check which suffixes still need icg_id
        need_lookup = [r for r in model_rows
                       if args.force or not r.get('icg_id')]
        if not need_lookup:
            # All rows already have icg_id; still refresh meta if forced
            if args.force:
                pass   # fall through to detail fetch below
            else:
                continue

        # ── 1. Search ICG for this modeldb_id ──────────────────────────────────
        print(f'[{idx:3d}/{len(modeldb_ids)}] modeldb_id={mid}: searching ICG ... ',
              end='', flush=True)
        icg_channels = icg_search(mid, delay)

        if not icg_channels:
            suffixes = sorted({r['suffix'] for r in model_rows})
            print(f'no ICG results  ({", ".join(suffixes)})')
            unmatched_total += len(model_rows)
            continue

        print(f'{len(icg_channels)} ICG channel(s) found')

        # ── 2. Match each suffix to an ICG channel ────────────────────────────
        for r in model_rows:
            if not args.force and r.get('icg_id'):
                continue
            suf   = r.get('suffix', '')
            match = _match(suf, icg_channels)
            if match:
                r['icg_id'] = str(match['icg_id'])
                matched_total += 1
            else:
                r['icg_id'] = ''
                unmatched_total += 1
                print(f'    [WARN] no ICG match for suffix={suf!r}  '
                      f'(available: {[c["icg_suffix"] for c in icg_channels]})')

        # ── 3. Fetch detail for each newly matched icg_id ────────────────────
        seen_icg_ids = set()
        for r in model_rows:
            raw = r.get('icg_id', '')
            if not raw:
                continue
            try:
                icg_id = int(raw)
            except ValueError:
                continue

            if icg_id in seen_icg_ids:
                continue
            seen_icg_ids.add(icg_id)

            if not args.force and icg_id in meta:
                continue   # already cached

            ref = icg_detail(icg_id, delay)
            new_detail_total += 1

            # Build meta row
            meta[icg_id] = {
                'icg_id':          str(icg_id),
                'modeldb_id':        str(mid),
                'suffix':          r.get('suffix', ''),
                'ion_class':       r.get('ion_class', ''),
                'icg_ref_id':      str(ref.get('id', '')),
                'authors':         ref.get('authors', ''),
                'title':           ref.get('title', ''),
                'year':            str(ref.get('date', '')),
                'pubmedid':        str(ref.get('id_pubmed', '')),
                'citation_count':  str(ref.get('citations', '')),
            }
            cit = ref.get('citations', '?')
            print(f'    icg_id={icg_id}  suffix={r["suffix"]!r}'
                  f'  citations={cit}')

    print()
    print(f'Matched:   {matched_total}')
    print(f'Unmatched: {unmatched_total}')
    print(f'New detail fetches: {new_detail_total}')

    if not args.dry_run:
        write_channel_db(chan_db_path, fieldnames, rows)
        print(f'Wrote channel_db → {chan_db_path}')
        write_meta_csv(meta_csv_path, meta)
        print(f'Wrote icg_channel_meta → {meta_csv_path}')
    else:
        print('[dry-run] No files written.')


if __name__ == '__main__':
    main()
