#!/usr/bin/env python3
"""
update_citations_icg.py
=======================
Refresh citation counts (and ref metadata) in icg_channel_meta.csv by
fetching each channel's detail record from the ICGenealogy API using its
ICG ID directly.

Requires that channel_db.csv already has the `icg_id` column populated
(run build_icg_meta.py first).

Usage
-----
    python update_citations_icg.py [--channel-db PATH] [--meta-csv PATH]
                                   [--rate FLOAT] [--dry-run]
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ICG_BASE   = 'https://icg.neurotheory.ox.ac.uk'
DETAIL_URL = ICG_BASE + '/api/app/chs/{icg_id}/'

SCRIPT_DIR       = Path(__file__).parent
DEFAULT_CHAN_DB   = SCRIPT_DIR / 'data' / 'channel_db.csv'
DEFAULT_META_CSV  = SCRIPT_DIR / 'data' / 'icg_channel_meta.csv'

META_FIELDS = [
    'icg_id', 'modeldb_id', 'suffix', 'ion_class',
    'icg_ref_id', 'authors', 'title', 'year', 'pubmedid', 'citation_count',
]


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
        print(f'  [WARN] {e}', file=sys.stderr)
        return None


def load_icg_ids(path: Path) -> list:
    """Return unique {icg_id, modeldb_id, suffix, ion_class} entries from channel_db.csv."""
    seen = set()
    result = []
    with open(path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            raw = row.get('icg_id', '').strip()
            if not raw:
                continue
            try:
                icg_id = int(raw)
            except ValueError:
                continue
            if icg_id in seen:
                continue
            seen.add(icg_id)
            result.append({
                'icg_id':    icg_id,
                'modeldb_id': row.get('modeldb_id', ''),
                'suffix':    row.get('suffix', ''),
                'ion_class': row.get('ion_class', ''),
            })
    return result


def load_meta(path: Path) -> dict:
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


def write_meta(path: Path, meta: dict) -> None:
    rows = [meta[k] for k in sorted(meta)]
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=META_FIELDS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--channel-db', default=str(DEFAULT_CHAN_DB))
    ap.add_argument('--meta-csv',   default=str(DEFAULT_META_CSV))
    ap.add_argument('--rate', type=float, default=1.0,
                    help='API requests per second (default 1.0)')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    delay         = 1.0 / max(args.rate, 0.1)
    chan_db_path  = Path(args.channel_db)
    meta_csv_path = Path(args.meta_csv)

    channels = load_icg_ids(chan_db_path)
    meta     = load_meta(meta_csv_path)

    if not channels:
        print('No icg_id values found in channel_db.csv.')
        print('Run build_icg_meta.py first to populate ICG IDs.')
        sys.exit(1)

    print(f'Updating {len(channels)} unique ICG channels  '
          f'(rate={args.rate:.1f} req/s)')

    updated = 0
    missing = 0

    for i, ch in enumerate(channels, 1):
        icg_id = ch['icg_id']
        print(f'[{i:4d}/{len(channels)}] icg_id={icg_id} ... ', end='', flush=True)

        data = _fetch(DETAIL_URL.format(icg_id=icg_id))
        time.sleep(delay)

        if not isinstance(data, dict):
            print('not found')
            missing += 1
            continue

        ref = data.get('ref') or {}
        old_cit = meta.get(icg_id, {}).get('citation_count', '')
        new_cit = str(ref.get('citations', ''))

        existing = meta.get(icg_id, {})
        meta[icg_id] = {
            'icg_id':         str(icg_id),
            'modeldb_id':     ch['modeldb_id'],
            'suffix':         ch['suffix'],
            'ion_class':      ch['ion_class'],
            'icg_ref_id':     str(ref.get('id', existing.get('icg_ref_id', ''))),
            'authors':        ref.get('authors') or existing.get('authors', ''),
            'title':          ref.get('title') or existing.get('title', ''),
            'year':           str(ref.get('date', existing.get('year', ''))),
            'pubmedid':       str(ref.get('id_pubmed', existing.get('pubmedid', ''))),
            'citation_count': new_cit or old_cit,
        }
        updated += 1
        print(f'citations: {old_cit!r} → {new_cit!r}')

    print(f'\nUpdated: {updated}  Not found: {missing}')

    if not args.dry_run:
        write_meta(meta_csv_path, meta)
        print(f'Wrote {len(meta)} rows → {meta_csv_path}')
    else:
        print('[dry-run] No files written.')


if __name__ == '__main__':
    main()
