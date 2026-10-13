"""Shared helpers for the GMO loci pipeline (stdlib + Biopython only)."""
import os, re, sys, json, time, csv, urllib.request, urllib.parse
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('GMO_ROOT', os.path.abspath(os.path.join(HERE, '..', '..')))  # repository root (see gmo_paths.py)
CFG_PATH = os.environ.get('GMO_CONFIG', os.path.join(HERE, '..', 'config.yaml'))


def load_config(path=CFG_PATH):
    cfg = yaml.safe_load(open(path))
    P = cfg['project_dir'] = os.path.expandvars(cfg['project_dir'])
    def ap(x): return x if os.path.isabs(x) else os.path.join(P, x)
    cfg['_paths'] = {
        'blast_stats': ap(cfg['inputs']['blast_stats']),
        'rep_seqs_dir': ap(cfg['inputs']['rep_seqs_dir']),
        'curated_dir': ap(cfg['inputs']['curated_dir']),
        'quarantine_tsv': ap(cfg['inputs']['quarantine_tsv']),
        'out': ap(cfg['out_dir']),
    }
    o = cfg['_paths']['out']
    for sub in ['registry', 'contigs', 'guides', 'logs']:
        cfg['_paths'][sub] = os.path.join(o, sub)
    return cfg


def round_dir(cfg, rnd):
    return os.path.join(cfg['_paths']['out'], f'round{rnd}')


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- FASTA
def read_fasta(path):
    recs, h, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if h is not None: recs.append((h, ''.join(buf)))
                h, buf = line[1:], []
            elif line:
                buf.append(line)
    if h is not None: recs.append((h, ''.join(buf)))
    return recs


def write_fasta(path, recs, width=80):
    with open(path, 'w') as f:
        for h, s in recs:
            f.write(f'>{h}\n')
            for i in range(0, len(s), width): f.write(s[i:i+width] + '\n')


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))[::-1]


def gc(s):
    s = s.upper(); n = sum(s.count(c) for c in 'ACGT')
    return 100.0 * (s.count('G') + s.count('C')) / n if n else float('nan')


# ---------------------------------------------------------------- TSV
def read_tsv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, cols=None):
    if not rows and not cols: open(path, 'w').close(); return
    cols = cols or list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


# ---------------------------------------------------------------- NCBI E-utilities
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
_last = [0.0]


def _api_key(cfg):
    return os.environ.get(cfg.get('ncbi_api_key_env', 'NCBI_API_KEY'), '')


def eutils(cfg, tool, params, tries=6, timeout=300, as_json=False):
    """Rate-limited GET with retries. tool = 'esearch.fcgi' | 'esummary.fcgi' | 'efetch.fcgi' | 'elink.fcgi'."""
    key = _api_key(cfg)
    gap = 0.11 if key else 0.36
    if key: params = dict(params, api_key=key)
    url = E + tool + '?' + urllib.parse.urlencode(params)
    err = None
    for t in range(tries):
        wait = gap - (time.time() - _last[0])
        if wait > 0: time.sleep(wait)
        _last[0] = time.time()
        try:
            data = urllib.request.urlopen(url, timeout=timeout).read().decode()
            if as_json:
                return json.loads(data)
            return data
        except Exception as e:
            err = e; time.sleep(2 + 3 * t)
    raise RuntimeError(f'eutils failed: {tool} {params}: {err}')


def wgs_prefix(acc):
    """JAHRGV010000067 -> JAHRGV01 ; PVXX01007568 -> PVXX01"""
    m = re.match(r'^([A-Z]{6}\d{2}|[A-Z]{4}\d{2})', acc)
    return m.group(1) if m else None


IS_WGS_CONTIG = re.compile(r'^[A-Z]{4,6}\d{8,11}(\.\d+)?$')
