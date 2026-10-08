#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mcr_africa_rebuild.py

Reproducible rebuild of the mcr-positive Escherichia coli / Shigella cohort from
Africa held in NCBI Pathogen Detection.

Design principle: retrieve the FULL global per-isolate metadata table for one
pinned Pathogen Detection build, then do mcr selection and country assignment
LOCALLY. No country filtering is applied at query time, because that is the
failure mode this rebuild exists to correct.

Stages
  1  Pin the build, download the per-isolate metadata TSV from the Pathogen
     Detection FTP metadata directory, record URL / size / SHA-256, and capture
     the AMRFinderPlus and AMRFinder database versions NCBI states for the build.
  2  Select mcr-positive isolates by tokenised whole-token regex match, with
     qualified calls (PARTIAL, PARTIAL_END_OF_CONTIG, HMM, INTERNAL_STOP, POINT)
     counted and ruled on explicitly.
  3  Assign country locally from a 54-state African Union alias table, with
     Unicode normalisation, word-boundary matching, documented fallbacks and a
     full residual list.
  4  Resolve retained isolates to sequence and link SRA runs through E-utilities.
  5  Reconcile against the earlier 118-isolate cohort.
  6  Write six outputs plus a provenance log.

Usage
  python mcr_africa_rebuild.py --outdir ./mcr_rebuild \
      --prior prior_cohort_118.tsv \
      [--build PDG000000004.6345] [--email you@example.org] [--api-key KEY] \
      [--skip-eutils] [--reuse-metadata PATH]

Dependencies: Python 3.8+ standard library only. requests is used when present
but is not required.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import io
import json
import os
import re
import ssl
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime, timezone

VERSION = "1.0"
UA = "mcr-africa-rebuild/%s (+pathogen-detection rebuild; contact: %s)"

FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Escherichia_coli_Shigella/"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TAXGROUP_SHORT = "Escherichia_coli_Shigella"

# Increase the CSV field cap: AMR_genotypes fields are long.
csv.field_size_limit(min(sys.maxsize, 2147483647))


# --------------------------------------------------------------------------
# logging
# --------------------------------------------------------------------------

class Log:
    """Collects the provenance log required by Step 6.5 and echoes to stderr."""

    def __init__(self):
        self.lines = []
        self.counts = OrderedDict()
        self.decisions = OrderedDict()

    def say(self, msg):
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        line = "[%s] %s" % (stamp, msg)
        self.lines.append(line)
        sys.stderr.write(line + "\n")
        sys.stderr.flush()

    def count(self, key, value):
        """Every filtering step records a count. Step with no count is a bug."""
        self.counts[key] = value
        self.say("COUNT  %-52s %s" % (key, value))

    def decision(self, key, value):
        self.decisions[key] = value
        self.say("DECIDE %-52s %s" % (key, value))

    def dump(self, path, extra_sections=None):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("mcr_africa_rebuild.py provenance log\n")
            fh.write("script version: %s\n" % VERSION)
            fh.write("python: %s\n" % sys.version.replace("\n", " "))
            fh.write("run started (UTC): %s\n\n" % RUN_START)
            fh.write("=" * 78 + "\nCOUNTS AT EVERY FILTERING STEP\n" + "=" * 78 + "\n")
            for k, v in self.counts.items():
                fh.write("%-56s %s\n" % (k, v))
            fh.write("\n" + "=" * 78 + "\nDISCRETIONARY DECISIONS\n" + "=" * 78 + "\n")
            for k, v in self.decisions.items():
                fh.write("%-32s %s\n" % (k, v))
            if extra_sections:
                for title, body in extra_sections:
                    fh.write("\n" + "=" * 78 + "\n" + title + "\n" + "=" * 78 + "\n")
                    fh.write(body.rstrip() + "\n")
            fh.write("\n" + "=" * 78 + "\nEVENT LOG\n" + "=" * 78 + "\n")
            fh.write("\n".join(self.lines) + "\n")


RUN_START = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
LOG = Log()


# --------------------------------------------------------------------------
# http
# --------------------------------------------------------------------------

def _opener(email):
    ctx = ssl.create_default_context()
    op = urllib.request.build_opener(urllib.request.HTTPSHandler(context=ctx))
    op.addheaders = [("User-Agent", UA % (VERSION, email or "unset"))]
    return op


def http_get(url, email=None, retries=4, timeout=120):
    last = None
    for attempt in range(retries):
        try:
            with _opener(email).open(url, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            wait = 2 ** attempt
            LOG.say("HTTP retry %d/%d after %ss for %s (%s)" % (attempt + 1, retries, wait, url, exc))
            time.sleep(wait)
    raise RuntimeError("GET failed after %d attempts: %s (%s)" % (retries, url, last))


def http_download(url, dest, email=None, chunk=1 << 20):
    """Stream to disk and return (bytes, sha256, md5)."""
    sha, md5, total = hashlib.sha256(), hashlib.md5(), 0
    with _opener(email).open(url, timeout=600) as resp, open(dest, "wb") as out:
        while True:
            buf = resp.read(chunk)
            if not buf:
                break
            out.write(buf)
            sha.update(buf)
            md5.update(buf)
            total += len(buf)
    return total, sha.hexdigest(), md5.hexdigest()


def list_dir(url, email=None):
    """Parse an NCBI FTP-over-HTTPS autoindex page into entry names."""
    body = http_get(url, email=email).decode("utf-8", "replace")
    names = re.findall(r'href="([^"?][^"]*)"', body)
    out = []
    for n in names:
        n = html.unescape(n)
        if n.startswith("/") or n.startswith("http") or n in ("../",):
            continue
        out.append(n)
    return out


# --------------------------------------------------------------------------
# STEP 1 - pin the build and fetch the metadata table
# --------------------------------------------------------------------------

BUILD_RE = re.compile(r"^(PDG\d{9}\.\d+)/?$")


def discover_builds(email):
    entries = list_dir(FTP_ROOT, email=email)
    builds = []
    for e in entries:
        m = BUILD_RE.match(e)
        if m:
            builds.append(m.group(1))
    builds.sort(key=lambda b: int(b.split(".")[1]))
    return builds


def step1(args):
    LOG.say("STEP 1 - pin build and fetch metadata")
    if args.offline:
        if not (args.build and args.reuse_metadata):
            raise SystemExit("--offline requires both --build and --reuse-metadata")
        LOG.say("OFFLINE mode: no FTP listing, checksumming the supplied local table only")
        local = args.reuse_metadata
        nbytes = os.path.getsize(local)
        sha, md5 = hashlib.sha256(), hashlib.md5()
        with open(local, "rb") as fh:
            for buf in iter(lambda: fh.read(1 << 20), b""):
                sha.update(buf)
                md5.update(buf)
        LOG.decision("1.2 build accession", args.build)
        LOG.count("1.3 metadata file size (bytes)", nbytes)
        return {"build": args.build,
                "build_url": urllib.parse.urljoin(FTP_ROOT, args.build + "/"),
                "metadata_file": os.path.basename(local),
                "metadata_url": "(offline run; file supplied locally)",
                "metadata_local": local, "metadata_bytes": nbytes,
                "metadata_sha256": sha.hexdigest(), "metadata_md5": md5.hexdigest(),
                "metadata_last_modified": "(offline run)",
                "build_dir_entries": ["(offline run; not listed)"],
                "metadata_dir_entries": ["(offline run; not listed)"]}
    if args.build:
        build = args.build
        LOG.say("build pinned by --build: %s" % build)
    else:
        builds = discover_builds(args.email)
        if not builds:
            raise SystemExit("no PDG builds found under %s" % FTP_ROOT)
        build = builds[-1]
        LOG.say("discovered %d builds; latest is %s" % (len(builds), build))
    LOG.decision("1.2 build accession", build)

    build_url = urllib.parse.urljoin(FTP_ROOT, build + "/")
    top = list_dir(build_url, email=args.email)
    LOG.say("build directory entries: %s" % ", ".join(sorted(top)))

    meta_url = urllib.parse.urljoin(build_url, "Metadata/")
    meta_entries = list_dir(meta_url, email=args.email)
    LOG.say("Metadata/ entries: %s" % ", ".join(sorted(meta_entries)))

    # Build date: the Last-Modified of the metadata directory listing is not
    # reliable, so take the build's own stated date from the directory listing
    # table when present, else the metadata file's Last-Modified header.
    candidates = [e for e in meta_entries if e.endswith(".metadata.tsv")]
    if not candidates:
        candidates = [e for e in meta_entries if e.endswith(".metadata.tsv.gz")]
    if not candidates:
        raise SystemExit("no *.metadata.tsv under %s (saw: %s)" % (meta_url, meta_entries))
    # prefer the file named for this build
    candidates.sort(key=lambda e: (not e.startswith(build), len(e)))
    fname = candidates[0]
    file_url = urllib.parse.urljoin(meta_url, fname)

    os.makedirs(args.outdir, exist_ok=True)
    local = args.reuse_metadata or os.path.join(args.outdir, fname)

    if args.reuse_metadata:
        LOG.say("reusing local metadata file: %s" % local)
        nbytes = os.path.getsize(local)
        sha, md5 = hashlib.sha256(), hashlib.md5()
        with open(local, "rb") as fh:
            for buf in iter(lambda: fh.read(1 << 20), b""):
                sha.update(buf)
                md5.update(buf)
        sha, md5 = sha.hexdigest(), md5.hexdigest()
        lastmod = "(reused local file; not re-fetched)"
    else:
        LOG.say("downloading %s" % file_url)
        try:
            with _opener(args.email).open(file_url, timeout=120) as r:
                lastmod = r.headers.get("Last-Modified", "(not reported)")
        except Exception:  # noqa: BLE001
            lastmod = "(not reported)"
        nbytes, sha, md5 = http_download(file_url, local, email=args.email)

    LOG.say("metadata file: %s" % fname)
    LOG.say("metadata url: %s" % file_url)
    LOG.say("metadata bytes: %d" % nbytes)
    LOG.say("metadata sha256: %s" % sha)
    LOG.say("metadata md5: %s" % md5)
    LOG.say("metadata Last-Modified: %s" % lastmod)
    LOG.count("1.3 metadata file size (bytes)", nbytes)

    return {
        "build": build,
        "build_url": build_url,
        "metadata_file": fname,
        "metadata_url": file_url,
        "metadata_local": local,
        "metadata_bytes": nbytes,
        "metadata_sha256": sha,
        "metadata_md5": md5,
        "metadata_last_modified": lastmod,
        "build_dir_entries": sorted(top),
        "metadata_dir_entries": sorted(meta_entries),
    }


def open_table(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def pick_column(header, candidates, required=True, what=""):
    lower = {h.lower(): h for h in header}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise SystemExit("could not find a column for %s; tried %s; header has %d columns: %s"
                         % (what, candidates, len(header), header[:40]))
    return None


# --------------------------------------------------------------------------
# STEP 2 - mcr selection
# --------------------------------------------------------------------------

MCR_RE = re.compile(r"^mcr-[0-9]+(\.[0-9]+)?$", re.IGNORECASE)

# AMRFinderPlus qualifiers seen appended to a gene symbol in the Pathogen
# Detection AMR_genotypes field, normally as SYMBOL=QUALIFIER.
QUALIFIERS = ("PARTIAL_END_OF_CONTIG", "PARTIAL_CONTIG_END", "PARTIAL",
              "INTERNAL_STOP", "HMM", "POINT", "MISTRANSLATION", "FRAME_SHIFT")

# Step 2.3 ruling. Changeable with --keep / --drop but these are the defaults
# used for the reported cohort.
DEFAULT_KEEP = {"COMPLETE", "PARTIAL", "PARTIAL_END_OF_CONTIG", "PARTIAL_CONTIG_END", "INTERNAL_STOP"}
DEFAULT_DROP = {"HMM", "POINT", "MISTRANSLATION", "FRAME_SHIFT"}

STEP23_RATIONALE = (
    "COMPLETE (no qualifier) kept: full-length nucleotide-identity call.\n"
    "PARTIAL and PARTIAL_END_OF_CONTIG kept: a gene truncated by a contig break "
    "is still evidence that the assembly carries mcr, which is what a carriage "
    "survey counts. Draft assemblies from short reads routinely break plasmid-borne "
    "mcr across contigs, so dropping these would bias the cohort against the very "
    "genomes most likely to carry the gene on a mobile element. Every such isolate "
    "is flagged in the output so the decision can be reversed downstream.\n"
    "INTERNAL_STOP kept but flagged: the gene is present; whether the allele is "
    "functional is a separate question that sequence inspection, not the metadata "
    "table, should answer.\n"
    "HMM excluded: an HMM-only hit has no allele-level nucleotide support and "
    "cannot be assigned an mcr variant, so it fails the allele regex in spirit as "
    "well as in letter.\n"
    "POINT excluded: POINT entries are point mutations in core genes; mcr is an "
    "acquired gene and a POINT-qualified mcr call would be a schema anomaly. The "
    "count is reported so any such anomaly is visible rather than silent."
)


NULL_TOKENS = {"NULL", "NA", "N/A", "NONE", "MISSING", "NOT COLLECTED", "NOT APPLICABLE",
               "NOT PROVIDED", "NOT DETERMINED", "UNKNOWN", "UNCALCULATED", "-", ""}


def clean(v):
    """NCBI writes several spellings of 'no value' into these columns. Treat them
    all as empty so they never become a deduplication key, a run accession or a
    sector signal."""
    s = (v or "").strip()
    return "" if s.upper() in NULL_TOKENS else s


def split_genotype(field):
    """Step 2.2 - treat the AMR genotype field as a delimited list."""
    if not field:
        return []
    parts = re.split(r"[,;]", field)
    return [p.strip() for p in parts if p.strip()]


def parse_token(token):
    """Return (base_symbol, qualifier). Qualifier is COMPLETE when absent."""
    tok = token.strip()
    qual = "COMPLETE"
    if "=" in tok:
        base, _, rest = tok.partition("=")
        rest_u = rest.strip().upper()
        for q in QUALIFIERS:
            if rest_u == q or rest_u.startswith(q):
                qual = q
                break
        else:
            qual = rest_u or "COMPLETE"
        tok = base.strip()
    else:
        up = tok.upper()
        for q in QUALIFIERS:
            if up.endswith("_" + q):
                qual = q
                tok = tok[: -(len(q) + 1)]
                break
    return tok.strip(), qual


def step2(meta, args):
    LOG.say("STEP 2 - select mcr-positive isolates")
    keep = set(args.keep) if args.keep else set(DEFAULT_KEEP)
    drop = set(args.drop) if args.drop else set(DEFAULT_DROP)
    LOG.decision("2.3 qualifiers kept", ", ".join(sorted(keep)))
    LOG.decision("2.3 qualifiers excluded", ", ".join(sorted(drop)))

    rows_out = []
    qual_counter = Counter()          # qualifier -> count of mcr TOKENS
    qual_isolates = defaultdict(set)  # qualifier -> set of isolate keys
    variant_counter = Counter()
    total = 0
    header = None
    cols = {}

    with open_table(meta["metadata_local"]) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        header = rdr.fieldnames or []
        LOG.say("metadata columns (%d): %s" % (len(header), ", ".join(header)))

        cols["amr"] = pick_column(header, ["AMR_genotypes", "amr_genotypes"], what="AMR genotype")
        cols["amr_core"] = pick_column(header, ["AMR_genotypes_core"], required=False)
        cols["target"] = pick_column(header, ["target_acc"], required=False)
        cols["biosample"] = pick_column(header, ["biosample_acc", "biosample"], what="BioSample")
        cols["asm"] = pick_column(header, ["asm_acc", "assembly_acc"], required=False)
        cols["run"] = pick_column(header, ["Run", "run_acc", "sra_run"], required=False)
        cols["platform"] = pick_column(header, ["Platform"], required=False)
        cols["layout"] = pick_column(header, ["LibraryLayout"], required=False)
        cols["geo"] = pick_column(header, ["geo_loc_name", "geo_loc_name_country", "location"], what="location")
        cols["latlon"] = pick_column(header, ["lat_lon", "latitude_longitude"], required=False)
        cols["host"] = pick_column(header, ["host"], required=False)
        cols["source"] = pick_column(header, ["isolation_source"], required=False)
        cols["coll"] = pick_column(header, ["collection_date"], required=False)
        cols["strain"] = pick_column(header, ["strain", "isolate"], required=False)
        cols["org"] = pick_column(header, ["scientific_name", "organism"], required=False)
        cols["bioproject"] = pick_column(header, ["bioproject_acc"], required=False)
        cols["created"] = pick_column(header, ["target_creation_date", "creation_date"], required=False)
        cols["epi"] = pick_column(header, ["epi_type"], required=False)
        cols["taxgroup"] = pick_column(header, ["taxgroup_name"], required=False)

        # Step 1.4 - versions NCBI states, harvested from the table itself.
        version_cols = [h for h in header if re.search(r"(amrfinder|refgene|db_version|_version)", h, re.I)]
        version_values = defaultdict(Counter)
        LOG.say("version-bearing columns: %s" % (", ".join(version_cols) or "(none)"))

        for row in rdr:
            total += 1
            for vc in version_cols:
                v = (row.get(vc) or "").strip()
                if v and v.upper() not in ("NULL", "NA", ""):
                    version_values[vc][v] += 1

            field = row.get(cols["amr"]) or ""
            if "mcr" not in field.lower():
                continue

            hits = []
            for tok in split_genotype(field):
                base, qual = parse_token(tok)
                if MCR_RE.match(base):
                    hits.append((base, qual))
            if not hits:
                continue

            kept_hits = [(b, q) for (b, q) in hits if q in keep]
            key = row.get(cols["biosample"]) or row.get(cols["target"]) or ""
            for b, q in hits:
                qual_counter[q] += 1
                qual_isolates[q].add(key)
            if not kept_hits:
                continue
            for b, _q in kept_hits:
                variant_counter[b.lower()] += 1

            def g(key):
                c = cols.get(key)
                return clean(row.get(c)) if c else ""

            rec = {
                "pd_target": g("target"),
                "biosample": g("biosample"),
                "assembly": g("asm"),
                "run_metadata": g("run"),
                "platform_metadata": g("platform"),
                "layout_metadata": g("layout"),
                "organism": g("org"),
                "taxgroup": g("taxgroup"),
                "strain": g("strain"),
                "host": g("host"),
                "isolation_source": g("source"),
                "collection_date": g("coll"),
                "location_string": g("geo"),
                "lat_lon": g("latlon"),
                "bioproject": g("bioproject"),
                "submission_date": g("created"),
                "epi_type": g("epi"),
                "amr_genotype": field,
                "mcr_variants": ";".join(sorted({b.lower() for b, _ in kept_hits})),
                "mcr_qualifiers": ";".join("%s=%s" % (b.lower(), q) for b, q in sorted(kept_hits)),
                "mcr_has_partial": "yes" if any(q.startswith("PARTIAL") for _, q in kept_hits) else "no",
                "mcr_has_internal_stop": "yes" if any(q == "INTERNAL_STOP" for _, q in kept_hits) else "no",
                "mcr_excluded_calls": ";".join("%s=%s" % (b.lower(), q) for b, q in hits if q not in keep) or "",
            }
            rows_out.append(rec)

    LOG.count("2.0 rows in metadata table", total)
    for q in sorted(set(list(qual_counter.keys()) + list(keep) + list(drop))):
        LOG.count("2.3 mcr token qualifier %-26s tokens" % q, qual_counter.get(q, 0))
        LOG.count("2.3 mcr token qualifier %-26s isolates" % q, len(qual_isolates.get(q, ())))
    LOG.count("2.4 GLOBAL mcr-positive isolates (pre-geography)", len(rows_out))
    for v, n in sorted(variant_counter.items()):
        LOG.count("2.4 global isolates carrying %s" % v, n)

    vsummary = {}
    for vc, ctr in version_values.items():
        vsummary[vc] = ctr.most_common(10)
        LOG.say("1.4 %s: %s" % (vc, "; ".join("%s (n=%d)" % (k, n) for k, n in ctr.most_common(5))))
    if not vsummary:
        LOG.say("1.4 WARNING no AMRFinderPlus/database version column found in the metadata table. "
                "Record the versions manually from the build's AMRFinderPlus directory listing in this log.")

    return rows_out, cols, header, vsummary


# --------------------------------------------------------------------------
# STEP 3 - local country assignment
# --------------------------------------------------------------------------

def norm(s):
    """NFKD, strip accents, lowercase, punctuation to space, collapse space."""
    if s is None:
        return ""
    s = unicodedata.normalize("NFKD", str(s))
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    s = s.replace("’", "'").replace("‘", "'")
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


# 54 African Union member states. Each entry: canonical -> (alpha2, alpha3, [aliases])
# Aliases include official names, short names, historical names and spelling and
# punctuation variants. Normalisation is applied to both sides at match time, so
# accents and punctuation here are cosmetic.
AU_STATES = OrderedDict([
    ("Algeria", ("DZ", "DZA", ["Algeria", "People's Democratic Republic of Algeria", "Algerie", "Algérie"])),
    ("Angola", ("AO", "AGO", ["Angola", "Republic of Angola"])),
    ("Benin", ("BJ", "BEN", ["Benin", "Republic of Benin", "Dahomey"])),
    ("Botswana", ("BW", "BWA", ["Botswana", "Republic of Botswana", "Bechuanaland"])),
    ("Burkina Faso", ("BF", "BFA", ["Burkina Faso", "Burkina", "Upper Volta", "Haute Volta"])),
    ("Burundi", ("BI", "BDI", ["Burundi", "Republic of Burundi"])),
    ("Cabo Verde", ("CV", "CPV", ["Cabo Verde", "Cape Verde", "Republic of Cabo Verde", "Republic of Cape Verde"])),
    ("Cameroon", ("CM", "CMR", ["Cameroon", "Republic of Cameroon", "Cameroun"])),
    ("Central African Republic", ("CF", "CAF", [
        "Central African Republic", "CAR", "Centrafrique", "Republique Centrafricaine",
        "République Centrafricaine", "Central African Rep", "Central African Rep."])),
    ("Chad", ("TD", "TCD", ["Chad", "Republic of Chad", "Tchad"])),
    ("Comoros", ("KM", "COM", ["Comoros", "The Comoros", "Union of the Comoros", "Comores"])),
    ("Democratic Republic of the Congo", ("CD", "COD", [
        "Democratic Republic of the Congo", "Democratic Republic of Congo", "DR Congo", "DRC",
        "D.R. Congo", "DR_Congo", "Congo Dem Rep", "Congo, Dem. Rep.", "Congo Democratic Republic",
        "Congo-Kinshasa", "Congo Kinshasa", "Zaire", "Zaïre", "Republique Democratique du Congo"])),
    ("Republic of the Congo", ("CG", "COG", [
        "Republic of the Congo", "Republic of Congo", "Congo-Brazzaville", "Congo Brazzaville",
        "Congo, Rep.", "Congo Rep", "Congo"])),
    ("Cote d'Ivoire", ("CI", "CIV", [
        "Cote d'Ivoire", "Côte d'Ivoire", "Cote dIvoire", "Cote d Ivoire", "Ivory Coast",
        "Republic of Cote d'Ivoire", "CIV"])),
    ("Djibouti", ("DJ", "DJI", ["Djibouti", "Republic of Djibouti"])),
    ("Egypt", ("EG", "EGY", ["Egypt", "Arab Republic of Egypt", "Misr", "United Arab Republic"])),
    ("Equatorial Guinea", ("GQ", "GNQ", ["Equatorial Guinea", "Republic of Equatorial Guinea", "Guinea Ecuatorial", "Spanish Guinea"])),
    ("Eritrea", ("ER", "ERI", ["Eritrea", "State of Eritrea"])),
    ("Eswatini", ("SZ", "SWZ", ["Eswatini", "Kingdom of Eswatini", "Swaziland"])),
    ("Ethiopia", ("ET", "ETH", ["Ethiopia", "Federal Democratic Republic of Ethiopia", "Abyssinia"])),
    ("Gabon", ("GA", "GAB", ["Gabon", "Gabonese Republic"])),
    ("Gambia", ("GM", "GMB", ["Gambia", "The Gambia", "Republic of the Gambia", "Gambia, The"])),
    ("Ghana", ("GH", "GHA", ["Ghana", "Republic of Ghana", "Gold Coast"])),
    ("Guinea", ("GN", "GIN", ["Guinea", "Republic of Guinea", "Guinea-Conakry", "Guinea Conakry", "French Guinea"])),
    ("Guinea-Bissau", ("GW", "GNB", [
        "Guinea-Bissau", "Guinea Bissau", "GuineaBissau", "Republic of Guinea-Bissau", "Portuguese Guinea"])),
    ("Kenya", ("KE", "KEN", ["Kenya", "Republic of Kenya"])),
    ("Lesotho", ("LS", "LSO", ["Lesotho", "Kingdom of Lesotho", "Basutoland"])),
    ("Liberia", ("LR", "LBR", ["Liberia", "Republic of Liberia"])),
    ("Libya", ("LY", "LBY", [
        "Libya", "State of Libya", "Libyan Arab Jamahiriya", "Socialist People's Libyan Arab Jamahiriya", "Libyan Arab Republic"])),
    ("Madagascar", ("MG", "MDG", ["Madagascar", "Republic of Madagascar", "Malagasy Republic"])),
    ("Malawi", ("MW", "MWI", ["Malawi", "Republic of Malawi", "Nyasaland"])),
    ("Mali", ("ML", "MLI", ["Mali", "Republic of Mali", "French Sudan"])),
    ("Mauritania", ("MR", "MRT", ["Mauritania", "Islamic Republic of Mauritania", "Mauritanie"])),
    ("Mauritius", ("MU", "MUS", ["Mauritius", "Republic of Mauritius", "Maurice"])),
    ("Morocco", ("MA", "MAR", ["Morocco", "Kingdom of Morocco", "Maroc"])),
    ("Mozambique", ("MZ", "MOZ", ["Mozambique", "Republic of Mozambique", "Mocambique", "Moçambique", "Portuguese East Africa"])),
    ("Namibia", ("NA", "NAM", ["Namibia", "Republic of Namibia", "South West Africa", "South-West Africa"])),
    ("Niger", ("NE", "NER", ["Niger", "Republic of the Niger", "Republic of Niger"])),
    ("Nigeria", ("NG", "NGA", ["Nigeria", "Federal Republic of Nigeria", "Niger Delta"])),
    ("Rwanda", ("RW", "RWA", ["Rwanda", "Republic of Rwanda", "Ruanda"])),
    ("Sahrawi Arab Democratic Republic", ("EH", "ESH", [
        "Western Sahara", "Sahrawi Arab Democratic Republic", "SADR", "Sahrawi Republic", "Saharawi", "Rio de Oro"])),
    ("Sao Tome and Principe", ("ST", "STP", [
        "Sao Tome and Principe", "São Tomé and Príncipe", "Sao Tome & Principe", "Sao Tome",
        "Democratic Republic of Sao Tome and Principe"])),
    ("Senegal", ("SN", "SEN", ["Senegal", "Republic of Senegal", "Sénégal"])),
    ("Seychelles", ("SC", "SYC", ["Seychelles", "Republic of Seychelles"])),
    ("Sierra Leone", ("SL", "SLE", ["Sierra Leone", "Republic of Sierra Leone"])),
    ("Somalia", ("SO", "SOM", ["Somalia", "Federal Republic of Somalia", "Somaliland", "Somali Republic", "Soomaaliya"])),
    ("South Africa", ("ZA", "ZAF", ["South Africa", "Republic of South Africa", "RSA", "Union of South Africa"])),
    ("South Sudan", ("SS", "SSD", ["South Sudan", "Republic of South Sudan"])),
    ("Sudan", ("SD", "SDN", ["Sudan", "Republic of the Sudan", "Republic of Sudan", "Anglo-Egyptian Sudan"])),
    ("Tanzania", ("TZ", "TZA", [
        "Tanzania", "United Republic of Tanzania", "Zanzibar", "Tanganyika", "Tanzania, United Republic of"])),
    ("Togo", ("TG", "TGO", ["Togo", "Togolese Republic", "Togoland"])),
    ("Tunisia", ("TN", "TUN", ["Tunisia", "Republic of Tunisia", "Tunisie"])),
    ("Uganda", ("UG", "UGA", ["Uganda", "Republic of Uganda"])),
    ("Zambia", ("ZM", "ZMB", ["Zambia", "Republic of Zambia", "Northern Rhodesia"])),
    ("Zimbabwe", ("ZW", "ZWE", ["Zimbabwe", "Republic of Zimbabwe", "Rhodesia", "Southern Rhodesia"])),
])

# Territories that sit near or in the region but are not African Union member
# states. Matching one of these excludes the isolate, and the removal is counted.
NON_AU_EXCLUSIONS = OrderedDict([
    ("Reunion", ["Reunion", "Réunion", "La Reunion"]),
    ("Mayotte", ["Mayotte"]),
    ("Saint Helena", ["Saint Helena", "St Helena", "St. Helena"]),
    ("Ascension", ["Ascension Island", "Ascension"]),
    ("Tristan da Cunha", ["Tristan da Cunha"]),
    ("Canary Islands", ["Canary Islands", "Islas Canarias", "Canarias", "Tenerife", "Gran Canaria"]),
    ("Ceuta", ["Ceuta"]),
    ("Melilla", ["Melilla"]),
    ("Madeira", ["Madeira", "Funchal"]),
])

# Free-text phrases that are regional, not national, and must never resolve to a
# country. "southern Africa" in particular must not become South Africa.
# Phrases outside Africa, or not countries at all, that contain an African
# country name as a whole token. Word-boundary matching alone does not stop
# "Papua New Guinea" resolving to Guinea, so these are tested first and block
# the match. Each entry is checked against the whole candidate string.
NON_AFRICAN_GUARDS = [
    "papua new guinea", "new guinea", "guinea pig", "guinea fowl", "gulf of guinea",
    "french guiana", "guyana", "new caledonia", "georgia usa",
]

NON_SPECIFIC = [
    "africa", "south africa region", "southern africa", "sub saharan africa", "subsaharan africa",
    "west africa", "western africa", "east africa", "eastern africa", "north africa", "northern africa",
    "central africa", "sahel", "horn of africa", "not applicable", "not collected", "missing",
    "unknown", "none", "null", "na", "n a", "not provided", "restricted access", "uncalculated",
]

# ISO codes are matched only as standalone upper-case tokens, because two-letter
# codes collide with ordinary words once lowercased.
_ALPHA2 = {v[0]: k for k, v in AU_STATES.items()}
_ALPHA3 = {v[1]: k for k, v in AU_STATES.items()}


def build_alias_index():
    """Longest-first alias list so Equatorial Guinea beats Guinea."""
    idx = []
    for canon, (a2, a3, aliases) in AU_STATES.items():
        seen = set()
        for al in aliases:
            n = norm(al)
            if n and n not in seen:
                seen.add(n)
                idx.append((n, canon, al))
    idx.sort(key=lambda t: -len(t[0]))
    return idx


ALIAS_INDEX = build_alias_index()

EXCL_INDEX = []
for terr, names in NON_AU_EXCLUSIONS.items():
    for nme in names:
        EXCL_INDEX.append((norm(nme), terr, nme))
EXCL_INDEX.sort(key=lambda t: -len(t[0]))


def phrase_in(haystack_norm, needle_norm):
    """Step 3.4 - whole-token match with word boundaries on normalised text."""
    if not needle_norm:
        return False
    return re.search(r"(?:^|\s)" + re.escape(needle_norm) + r"(?:\s|$)", haystack_norm) is not None


def candidate_strings(location):
    """Step 3.2 - the whole string, plus the 'Country: subregion' convention,
    plus the reverse order some records use."""
    out = []
    s = (location or "").strip()
    if not s:
        return out
    out.append(s)
    if ":" in s:
        head, _, tail = s.partition(":")
        out.append(head)
        out.append(tail)
    for sep in (",", "/", "|"):
        if sep in s:
            parts = [p for p in s.split(sep) if p.strip()]
            out.extend(parts)
            if len(parts) >= 2:
                out.append(parts[-1])
    return out


def match_country(location, alias_hits=None, excl_hits=None):
    """Return (canonical_country or None, matched_alias or None, reason)."""
    if not location or not location.strip():
        return None, None, "empty location field"
    whole = norm(location)
    if whole in NON_SPECIFIC:
        return None, None, "non-specific location: %r" % location

    cands = candidate_strings(location)
    cand_norms = [norm(c) for c in cands if norm(c)]

    # Exclusions first: a Reunion record must not fall through to a country match.
    for n, terr, raw in EXCL_INDEX:
        for cn in cand_norms:
            if phrase_in(cn, n):
                if excl_hits is not None:
                    excl_hits[terr] += 1
                return None, raw, "excluded non-AU territory: %s" % terr

    # "southern africa" must be rejected before "south africa" can match inside it.
    for cn in cand_norms:
        if phrase_in(cn, "southern africa"):
            return None, None, "non-specific region: southern Africa"

    # Guard phrases that contain an African country name but are not that country.
    for g in NON_AFRICAN_GUARDS:
        if phrase_in(whole, g):
            return None, None, "non-African guard phrase matched: %s" % g

    for n, canon, raw in ALIAS_INDEX:
        for cn in cand_norms:
            if phrase_in(cn, n):
                if alias_hits is not None:
                    alias_hits[(canon, raw)] += 1
                return canon, raw, "alias match"

    # ISO alpha-3 is accepted ONLY when a candidate string is exactly the code,
    # never as a token inside a longer string: "Viet Nam" contains the token
    # "nam" and would otherwise resolve to Namibia. ISO alpha-2 is not matched at
    # all, because two-letter codes collide with US state abbreviations, so
    # "USA:GA" would resolve to Gabon, "USA:SC" to Seychelles and so on.
    for cn in cand_norms:
        for code, canon in _ALPHA3.items():
            if cn == norm(code):
                if alias_hits is not None:
                    alias_hits[(canon, code)] += 1
                return canon, code, "ISO alpha-3 exact-string match"

    return None, None, "no alias matched"


# Coarse bounding boxes for the lat/long fallback (Step 3.5a). Deliberately
# coarse and always flagged: a box hit is a lead for manual confirmation, never a
# silent assignment. Boxes are (min_lat, max_lat, min_lon, max_lon).
BBOX = {
    "Algeria": (18.9, 37.1, -8.7, 12.0), "Angola": (-18.1, -4.4, 11.6, 24.1),
    "Benin": (6.2, 12.4, 0.7, 3.9), "Botswana": (-26.9, -17.8, 19.9, 29.4),
    "Burkina Faso": (9.4, 15.1, -5.6, 2.4), "Burundi": (-4.5, -2.3, 29.0, 30.9),
    "Cabo Verde": (14.8, 17.2, -25.4, -22.7), "Cameroon": (1.6, 13.1, 8.4, 16.2),
    "Central African Republic": (2.2, 11.0, 14.4, 27.5), "Chad": (7.4, 23.5, 13.5, 24.0),
    "Comoros": (-12.5, -11.3, 43.2, 44.6), "Democratic Republic of the Congo": (-13.5, 5.4, 12.2, 31.3),
    "Republic of the Congo": (-5.1, 3.7, 11.1, 18.7), "Cote d'Ivoire": (4.3, 10.8, -8.6, -2.5),
    "Djibouti": (10.9, 12.8, 41.7, 43.5), "Egypt": (22.0, 31.7, 24.7, 36.9),
    "Equatorial Guinea": (0.9, 3.8, 5.6, 11.4), "Eritrea": (12.3, 18.1, 36.4, 43.2),
    "Eswatini": (-27.4, -25.7, 30.8, 32.2), "Ethiopia": (3.4, 14.9, 32.9, 48.0),
    "Gabon": (-4.0, 2.4, 8.6, 14.6), "Gambia": (13.0, 13.9, -17.0, -13.7),
    "Ghana": (4.7, 11.2, -3.3, 1.2), "Guinea": (7.1, 12.7, -15.1, -7.6),
    "Guinea-Bissau": (10.8, 12.7, -16.8, -13.6), "Kenya": (-4.7, 5.5, 33.9, 41.9),
    "Lesotho": (-30.7, -28.5, 27.0, 29.5), "Liberia": (4.3, 8.6, -11.5, -7.3),
    "Libya": (19.5, 33.2, 9.3, 25.2), "Madagascar": (-25.7, -11.9, 43.2, 50.5),
    "Malawi": (-17.2, -9.3, 32.6, 35.9), "Mali": (10.1, 25.0, -12.3, 4.3),
    "Mauritania": (14.7, 27.3, -17.1, -4.8), "Mauritius": (-20.6, -19.9, 57.2, 57.9),
    "Morocco": (27.6, 35.9, -13.2, -1.0), "Mozambique": (-26.9, -10.5, 30.2, 40.9),
    "Namibia": (-28.9, -16.9, 11.7, 25.3), "Niger": (11.7, 23.5, 0.2, 16.0),
    "Nigeria": (4.3, 13.9, 2.7, 14.7), "Rwanda": (-2.8, -1.1, 28.9, 30.9),
    "Sahrawi Arab Democratic Republic": (20.8, 27.7, -17.1, -8.7),
    "Sao Tome and Principe": (-0.1, 1.8, 6.4, 7.5), "Senegal": (12.3, 16.7, -17.6, -11.4),
    "Seychelles": (-10.3, -3.7, 46.2, 56.3), "Sierra Leone": (6.9, 10.0, -13.4, -10.2),
    "Somalia": (-1.7, 12.0, 40.9, 51.4), "South Africa": (-34.9, -22.1, 16.4, 32.9),
    "South Sudan": (3.5, 12.3, 24.1, 35.9), "Sudan": (8.7, 22.2, 21.8, 38.6),
    "Tanzania": (-11.8, -0.9, 29.3, 40.5), "Togo": (6.1, 11.1, -0.2, 1.8),
    "Tunisia": (30.2, 37.6, 7.5, 11.6), "Uganda": (-1.5, 4.2, 29.5, 35.0),
    "Zambia": (-18.1, -8.2, 21.9, 33.7), "Zimbabwe": (-22.5, -15.6, 25.2, 33.1),
}


def parse_lat_lon(s):
    """Parse the BioSample lat_lon convention, e.g. '12.34 N 7.89 E'."""
    if not s:
        return None
    m = re.search(r"([-+]?\d+(?:\.\d+)?)\s*([NS])[,\s]+([-+]?\d+(?:\.\d+)?)\s*([EW])", s, re.I)
    if m:
        # The hemisphere letter carries the sign. Take the magnitude first so a
        # record written "-26.2 S" is not double-negated into the wrong hemisphere.
        lat = abs(float(m.group(1))) * (1 if m.group(2).upper() == "N" else -1)
        lon = abs(float(m.group(3))) * (1 if m.group(4).upper() == "E" else -1)
        return lat, lon
    m = re.search(r"^\s*([-+]?\d+(?:\.\d+)?)[,\s]+([-+]?\d+(?:\.\d+)?)\s*$", s)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None


_RG = None
_RG_TRIED = False


def _rg():
    """Optional offline reverse geocoder. Bounding boxes alone are too coarse to
    separate neighbours, so when the `reverse_geocoder` package is installed it is
    used first and the box check becomes a sanity check."""
    global _RG, _RG_TRIED
    if not _RG_TRIED:
        _RG_TRIED = True
        try:
            import reverse_geocoder  # noqa: F401
            _RG = reverse_geocoder
            LOG.say("3.5a reverse_geocoder is installed and will be used for lat/long resolution")
        except Exception:  # noqa: BLE001
            _RG = None
            LOG.say("3.5a reverse_geocoder not installed; lat/long falls back to coarse bounding boxes "
                    "and ambiguous points are left unresolved rather than guessed")
    return _RG


def country_from_latlon(s):
    pt = parse_lat_lon(s)
    if not pt:
        return None, None
    lat, lon = pt

    rg = _rg()
    if rg is not None:
        try:
            hit = rg.search((lat, lon), mode=1, verbose=False)[0]
            cc = (hit.get("cc") or "").upper()
            if cc in _ALPHA2:
                return _ALPHA2[cc], "reverse_geocoder cc=%s at %.5f,%.5f" % (cc, lat, lon)
            if cc:
                return None, "reverse_geocoder resolved to non-AU country %s at %.5f,%.5f" % (cc, lat, lon)
        except Exception as exc:  # noqa: BLE001
            LOG.say("3.5a reverse_geocoder failed for %.5f,%.5f: %s" % (lat, lon, exc))

    hits = [c for c, (a, b, d, e) in BBOX.items() if a <= lat <= b and d <= lon <= e]
    if len(hits) == 1:
        return hits[0], "%.5f,%.5f" % (lat, lon)
    if len(hits) > 1:
        return None, "ambiguous bbox: %s at %.5f,%.5f" % ("|".join(sorted(hits)), lat, lon)
    return None, "no African bbox contains %.5f,%.5f" % (lat, lon)


def self_test_boundaries():
    """Step 3.4 - these cases are asserted and the result is reported."""
    cases = [
        ("Niger", "Niger", True), ("Nigeria", "Nigeria", True),
        ("Niger: Niamey", "Niger", True), ("Nigeria: Lagos", "Nigeria", True),
        ("Guinea", "Guinea", True), ("Equatorial Guinea", "Equatorial Guinea", True),
        ("Guinea-Bissau", "Guinea-Bissau", True), ("Guinea Bissau", "Guinea-Bissau", True),
        ("Papua New Guinea", None, True), ("Chad", "Chad", True),
        ("sample from a Chadian household", None, True),
        ("Cote d'Ivoire", "Cote d'Ivoire", True), ("Côte d'Ivoire", "Cote d'Ivoire", True),
        ("Ivory Coast", "Cote d'Ivoire", True), ("CIV", "Cote d'Ivoire", True),
        ("Congo, Dem. Rep.", "Democratic Republic of the Congo", True),
        ("Congo-Kinshasa", "Democratic Republic of the Congo", True),
        ("Zaire", "Democratic Republic of the Congo", True),
        ("DRC", "Democratic Republic of the Congo", True),
        ("Congo-Brazzaville", "Republic of the Congo", True),
        ("Republic of the Congo", "Republic of the Congo", True),
        ("Central African Republic", "Central African Republic", True),
        ("Centrafrique", "Central African Republic", True),
        ("Swaziland", "Eswatini", True), ("Eswatini", "Eswatini", True),
        ("Cape Verde", "Cabo Verde", True), ("Cabo Verde", "Cabo Verde", True),
        ("United Republic of Tanzania", "Tanzania", True), ("Zanzibar", "Tanzania", True),
        ("The Gambia", "Gambia", True), ("Gambia", "Gambia", True),
        ("São Tomé and Príncipe", "Sao Tome and Principe", True),
        ("Western Sahara", "Sahrawi Arab Democratic Republic", True),
        ("South Africa: Gauteng", "South Africa", True),
        ("southern Africa", None, True), ("Southern Africa", None, True),
        ("Upper Volta", "Burkina Faso", True), ("Dahomey", "Benin", True),
        ("Nyasaland", "Malawi", True), ("Northern Rhodesia", "Zambia", True),
        ("Rhodesia", "Zimbabwe", True), ("Abyssinia", "Ethiopia", True),
        ("Somaliland", "Somalia", True), ("Libyan Arab Jamahiriya", "Libya", True),
        # ISO-code collisions found in the real build on 2026-10-07.
        ("Viet Nam", None, True), ("Viet Nam:Thai Binh", None, True),
        ("Viet Nam: Ho Chi Minh City", None, True),
        ("USA:GA", None, True), ("USA:SC", None, True), ("USA:SD", None, True),
        ("USA:NE", None, True), ("USA:TN", None, True), ("USA: Boston, MA", None, True),
        ("USA:MA,Boston", None, True),
        ("NAM", "Namibia", True), ("NGA", "Nigeria", True), ("ZAF", "South Africa", True),
        ("Papua New Guinea: Port Moresby", None, True), ("Gulf of Guinea", None, True),
        ("guinea fowl farm", None, True), ("Niger Delta", "Nigeria", True),
        ("Niger Delta, Nigeria", "Nigeria", True), ("26.2 S 28.0 E", None, True),
        ("Reunion", None, True), ("Mayotte", None, True), ("Canary Islands", None, True),
        ("Ceuta", None, True), ("Melilla", None, True), ("Madeira", None, True),
        ("Kinshasa: Democratic Republic of the Congo", "Democratic Republic of the Congo", True),
    ]
    results, failures = [], 0
    for loc, expected, _ in cases:
        got, alias, reason = match_country(loc)
        ok = (got == expected)
        if not ok:
            failures += 1
        results.append((loc, expected, got, alias, reason, "PASS" if ok else "FAIL"))
    return results, failures


def eutils_biosample_geo(biosamples, email, api_key, pause):
    """Step 3.5b - geographic location recorded on the linked BioSample."""
    out = {}
    if not biosamples:
        return out
    batch = 200
    accs = list(biosamples)
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        params = {"db": "biosample", "term": " OR ".join("%s[accn]" % a for a in chunk),
                  "retmax": str(len(chunk)), "retmode": "json", "tool": "mcr_africa_rebuild"}
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key
        try:
            js = json.loads(http_get(EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params), email=email))
            uids = js.get("esearchresult", {}).get("idlist", [])
        except Exception as exc:  # noqa: BLE001
            LOG.say("biosample esearch failed for a chunk: %s" % exc)
            uids = []
        time.sleep(pause)
        if not uids:
            continue
        p2 = {"db": "biosample", "id": ",".join(uids), "rettype": "full", "retmode": "xml",
              "tool": "mcr_africa_rebuild"}
        if email:
            p2["email"] = email
        if api_key:
            p2["api_key"] = api_key
        try:
            xml = http_get(EUTILS + "efetch.fcgi?" + urllib.parse.urlencode(p2), email=email).decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            LOG.say("biosample efetch failed for a chunk: %s" % exc)
            time.sleep(pause)
            continue
        for blk in re.findall(r"<BioSample\b.*?</BioSample>", xml, re.S):
            accm = re.search(r'accession="(SAM[A-Z0-9]+)"', blk)
            if not accm:
                continue
            geo = None
            for am in re.finditer(r'<Attribute[^>]*attribute_name="([^"]+)"[^>]*>(.*?)</Attribute>', blk, re.S):
                if am.group(1).strip().lower() in ("geo_loc_name", "geographic location", "country",
                                                   "geographic location (country and/or sea)"):
                    geo = html.unescape(am.group(2)).strip()
                    break
            if geo:
                out[accm.group(1)] = geo
        time.sleep(pause)
    return out


def step3(rows, args):
    LOG.say("STEP 3 - assign country locally")
    tests, failures = self_test_boundaries()
    LOG.count("3.4 boundary test cases run", len(tests))
    LOG.count("3.4 boundary test failures", failures)
    if failures:
        LOG.say("3.4 WARNING boundary tests failed; inspect the boundary test report before trusting counts")

    LOG.decision("3.7 Western Sahara", (
        "Included as an African Union member state under the name Sahrawi Arab Democratic Republic. "
        "The SADR has been an AU (formerly OAU) member since 1984, so a rebuild scoped to AU member "
        "states must include it. Isolates assigned there are flagged western_sahara_review=yes so the "
        "disputed status is visible to a reader, and they are never merged into Morocco."))

    alias_hits = Counter()
    excl_hits = Counter()
    resolved, residual = [], []
    reason_counter = Counter()

    for r in rows:
        country, alias, reason = match_country(r["location_string"], alias_hits, excl_hits)
        r["country_match_alias"] = alias or ""
        r["country_assignment_reason"] = reason
        r["location_fallback_flag"] = "no"
        if country:
            r["assigned_country"] = country
            resolved.append(r)
        else:
            r["assigned_country"] = ""
            reason_counter[reason.split(":")[0]] += 1
            residual.append(r)

    LOG.count("3.2 isolates assigned from the location field", len(resolved))
    LOG.count("3.2 isolates unresolved after the location field", len(residual))
    for terr in NON_AU_EXCLUSIONS:
        LOG.count("3.6 excluded by territory %-22s" % terr, excl_hits.get(terr, 0))
    LOG.count("3.6 excluded by non-AU territories, total", sum(excl_hits.values()))

    # 3.5a lat/long fallback
    latlon_assigned = 0
    still = []
    for r in residual:
        if r.get("lat_lon"):
            c, note = country_from_latlon(r["lat_lon"])
            if c:
                r["assigned_country"] = c
                r["location_fallback_flag"] = "yes-latlon"
                r["country_assignment_reason"] = "lat/long bounding box: %s" % note
                resolved.append(r)
                latlon_assigned += 1
                continue
            if note:
                r["country_assignment_reason"] += " | latlon: %s" % note
        still.append(r)
    residual = still
    LOG.count("3.5a isolates assigned by lat/long fallback", latlon_assigned)

    # 3.5b BioSample fallback
    biosample_assigned = 0
    if not args.skip_eutils and residual:
        want = [r["biosample"] for r in residual if r.get("biosample")]
        LOG.say("3.5b querying BioSample for %d unresolved isolates" % len(want))
        geo = eutils_biosample_geo(want, args.email, args.api_key, args.pause)
        still = []
        for r in residual:
            g = geo.get(r.get("biosample", ""))
            if g:
                c, alias, reason = match_country(g, alias_hits, excl_hits)
                if c:
                    r["assigned_country"] = c
                    r["location_fallback_flag"] = "yes-biosample"
                    r["country_match_alias"] = alias or ""
                    r["country_assignment_reason"] = "BioSample geo_loc_name %r -> %s" % (g, reason)
                    resolved.append(r)
                    biosample_assigned += 1
                    continue
                r["country_assignment_reason"] += " | biosample geo %r: %s" % (g, reason)
            still.append(r)
        residual = still
    else:
        LOG.say("3.5b skipped (--skip-eutils or nothing unresolved)")
    LOG.count("3.5b isolates assigned by BioSample fallback", biosample_assigned)

    african = [r for r in resolved if r["assigned_country"]]
    for r in african:
        r["western_sahara_review"] = "yes" if r["assigned_country"] == "Sahrawi Arab Democratic Republic" else "no"

    LOG.count("3.0 AFRICAN UNION isolates retained", len(african))
    LOG.count("3.8 unresolved-location residual isolates (worldwide)", len(residual))
    for k, v in reason_counter.most_common():
        LOG.count("3.x residual reason: %s" % k, v)

    return african, residual, alias_hits, excl_hits, tests


# --------------------------------------------------------------------------
# STEP 4 - resolve to sequence
# --------------------------------------------------------------------------

def eutils_sra_for_biosamples(biosamples, email, api_key, pause):
    """Step 4.2 - esearch + esummary, db=sra, on a single stated date."""
    runs = defaultdict(list)
    accs = [a for a in biosamples if a]
    batch = 100
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        params = {"db": "sra", "term": " OR ".join("%s[accn]" % a for a in chunk),
                  "retmax": "5000", "retmode": "json", "tool": "mcr_africa_rebuild"}
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key
        try:
            js = json.loads(http_get(EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params), email=email))
            uids = js.get("esearchresult", {}).get("idlist", [])
        except Exception as exc:  # noqa: BLE001
            LOG.say("sra esearch failed for a chunk: %s" % exc)
            uids = []
        time.sleep(pause)
        for j in range(0, len(uids), 200):
            sub = uids[j:j + 200]
            p2 = {"db": "sra", "id": ",".join(sub), "retmode": "json", "tool": "mcr_africa_rebuild"}
            if email:
                p2["email"] = email
            if api_key:
                p2["api_key"] = api_key
            try:
                js2 = json.loads(http_get(EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(p2), email=email))
            except Exception as exc:  # noqa: BLE001
                LOG.say("sra esummary failed for a chunk: %s" % exc)
                time.sleep(pause)
                continue
            res = js2.get("result", {})
            for uid in res.get("uids", []):
                rec = res.get(uid, {})
                expxml = rec.get("expxml", "") or ""
                runxml = rec.get("runs", "") or ""
                bs = re.search(r"<Biosample>(SAM[A-Z0-9]+)</Biosample>", expxml)
                if not bs:
                    bs = re.search(r"(SAM[ENД]?[A-Z]?\d+)", expxml)
                bsa = bs.group(1) if bs else ""
                plat = re.search(r'<Platform[^>]*>([^<]*)<', expxml)
                layout = "PAIRED" if "<PAIRED" in expxml else ("SINGLE" if "<SINGLE" in expxml else "")
                inst = re.search(r'instrument_model="([^"]*)"', expxml)
                for rm in re.finditer(r'<Run acc="([^"]+)"[^>]*total_spots="(\d*)"[^>]*total_bases="(\d*)"', runxml):
                    runs[bsa].append({
                        "run": rm.group(1),
                        "platform": (plat.group(1).strip() if plat else ""),
                        "instrument": (inst.group(1) if inst else ""),
                        "layout": layout,
                        "total_spots": rm.group(2),
                        "total_bases": rm.group(3),
                    })
            time.sleep(pause)
    return runs


VFDB_INVASION = ["ipaH", "ipaB", "ipaC", "ipaD", "icsA (virG)", "virB", "mxiA/mxiD", "spa15/spa32"]
VFDB_COMMENSAL = ["ompA", "csgA/csgB/csgD (csg operon)", "ecpA/ecpB/ecpC/ecpD/ecpE (ecp pilus)"]


def step4(african, args):
    LOG.say("STEP 4 - resolve isolates to sequence")
    query_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    LOG.decision("4.2 E-utilities query date (UTC)", query_date)

    if args.skip_eutils:
        LOG.say("4.2 DEVIATION: E-utilities not queried (--skip-eutils). Run accessions, platform and "
                "library layout are taken from the build's own Run, Platform and LibraryLayout columns, "
                "which are NCBI-computed for this build. Read counts (total_spots, total_bases) are not "
                "in the metadata table and are therefore MISSING; recover them with a live esummary run.")
        n_with = 0
        for r in african:
            r["sra_runs"] = r.get("run_metadata", "")
            if r["sra_runs"]:
                n_with += 1
            r["sra_platform"] = r.get("platform_metadata", "")
            r["sra_layout"] = r.get("layout_metadata", "")
            r["sra_total_spots"] = "NOT RETRIEVED (no E-utilities call)"
            r["sra_total_bases"] = "NOT RETRIEVED (no E-utilities call)"
            r["sra_source"] = "build metadata columns Run/Platform/LibraryLayout (no E-utilities call)"
        LOG.count("4.2 isolates with at least one run accession (metadata column)", n_with)
        LOG.count("4.2 isolates with no run accession", len(african) - n_with)
        LOG.decision("4.2 DEVIATION FROM SPEC", (
            "E-utilities esearch/esummary was not run; this environment has no network route to NCBI. "
            "Run accession, platform and library layout come from the pinned build's own columns. "
            "Read counts are absent. Step 4.2 is therefore satisfied in part only."))
    else:
        runmap = eutils_sra_for_biosamples([r["biosample"] for r in african],
                                           args.email, args.api_key, args.pause)
        n_with = 0
        for r in african:
            rs = runmap.get(r["biosample"], [])
            if rs:
                n_with += 1
            r["sra_runs"] = ";".join(x["run"] for x in rs)
            r["sra_platform"] = ";".join(sorted({x["platform"] for x in rs if x["platform"]}))
            r["sra_layout"] = ";".join(sorted({x["layout"] for x in rs if x["layout"]}))
            r["sra_total_spots"] = ";".join(x["total_spots"] for x in rs)
            r["sra_total_bases"] = ";".join(x["total_bases"] for x in rs)
            r["sra_source"] = "E-utilities esearch/esummary db=sra on %s" % query_date
        LOG.count("4.2 isolates with at least one linked SRA run", n_with)
        LOG.count("4.2 isolates with no linked SRA run", len(african) - n_with)

    # 4.3 independent deduplication, collisions reported not dropped
    collisions = []
    for field, label in (("biosample", "BioSample"), ("assembly", "assembly accession")):
        seen = defaultdict(list)
        for r in african:
            v = (r.get(field) or "").strip()
            if v:
                seen[v].append(r.get("pd_target") or r.get("biosample"))
        dupes = {k: v for k, v in seen.items() if len(v) > 1}
        LOG.count("4.3 distinct %s values" % label, len(seen))
        LOG.count("4.3 %s collisions" % label, len(dupes))
        for k, v in sorted(dupes.items()):
            collisions.append((label, k, ";".join(v)))

    runseen = defaultdict(list)
    for r in african:
        for run in (r.get("sra_runs") or "").replace(",", ";").split(";"):
            run = run.strip()
            if run:
                runseen[run].append(r.get("biosample"))
    rdupes = {k: v for k, v in runseen.items() if len(set(v)) > 1}
    LOG.count("4.3 distinct run accessions", len(runseen))
    LOG.count("4.3 run accession collisions", len(rdupes))
    for k, v in sorted(rdupes.items()):
        collisions.append(("run accession", k, ";".join(sorted(set(v)))))

    # 4.4 species confirmation downstream of the submitted name
    to_screen = []
    for r in african:
        name = (r.get("organism") or "").strip()
        n = norm(name)
        is_ecoli = n.startswith("escherichia coli") or n == "escherichia coli"
        r["submitted_species_is_ecoli"] = "yes" if is_ecoli else "no"
        if not is_ecoli:
            r["species_decision"] = "PENDING-VFDB-SCREEN"
            to_screen.append(r)
        else:
            r["species_decision"] = "retain (submitted name is E. coli; organism group membership consistent)"
    LOG.count("4.4 isolates submitted under a non-E. coli name", len(to_screen))
    LOG.decision("4.4 species confirmation", (
        "Isolates submitted as Shigella or any other non-E. coli name are not retained or dropped on the "
        "submitted name. They are written to the species-screen list with the VFDB invasion-locus panel "
        "(%s) and the commensal E. coli marker panel (%s). Each needs a blastn screen against VFDB on its "
        "own assembly, and a stated retain or exclude decision per isolate, before the cohort is final. "
        "This script does not fabricate that result: the decision column reads PENDING-VFDB-SCREEN until "
        "the screen is done." % (", ".join(VFDB_INVASION), ", ".join(VFDB_COMMENSAL))))

    # Sector, only where host or isolation source actually says something.
    for r in african:
        r["assigned_sector"] = derive_sector(r.get("host"), r.get("isolation_source"))
    sect = Counter(r["assigned_sector"] for r in african)
    for k, v in sect.most_common():
        LOG.count("4.x sector %s" % k, v)

    return collisions, to_screen, query_date


POULTRY_PAT = re.compile(r"\b(gallus|chicken|broiler|layer|poultry|hen|cockerel|chick|duck|turkey|"
                         r"guinea fowl|anas|meleagris)\b", re.I)
HUMAN_PAT = re.compile(r"\b(homo sapiens|human|patient)\b", re.I)
ENV_PAT = re.compile(r"\b(wastewater|waste water|sewage|river|lake|soil|bathing water|effluent|"
                     r"surface water|drain|abattoir environment|environmental swab)\b", re.I)
FOOD_PAT = re.compile(r"\b(meat|beef|pork|mutton|cheese|milk|sausage|food|carcass|retail|"
                      r"cucumber|lettuce|tomato|spinach|vegetable|salad|fruit|produce|herb)\b", re.I)
ANIMAL_PAT = re.compile(r"\b(bos taurus|cattle|cow|sus scrofa|pig|piglet|swine|ovis|goat|sheep|canis|dog|"
                        r"felis|cat|camel|horse|equus|starling|bird|rodent|rat|mouse|fish|"
                        r"papio|macaca|chlorocebus|gorilla|pan troglodytes|baboon|monkey|primate|"
                        r"bat|antelope|buffalo|donkey|rabbit)\b", re.I)


def derive_sector(host, source):
    """Never guess. Absent or uninformative metadata becomes Unspecified."""
    blob = " ".join(x for x in (host or "", source or "") if x and x.upper() not in ("NULL", "NA", "MISSING"))
    if not blob.strip():
        return "Unspecified"
    if HUMAN_PAT.search(blob):
        return "Human"
    if POULTRY_PAT.search(blob):
        return "Poultry"
    if ENV_PAT.search(blob):
        return "Environment"
    if ANIMAL_PAT.search(blob):
        return "Other animal"
    if FOOD_PAT.search(blob):
        return "Food (non-poultry)"
    return "Unspecified"


# --------------------------------------------------------------------------
# STEP 5 - reconcile against the earlier cohort
# --------------------------------------------------------------------------

KNOWN_MISSING_COUNTRIES = ["South Africa", "Central African Republic", "Democratic Republic of the Congo"]


def step5(african, prior_path, prior_query_date, args):
    LOG.say("STEP 5 - reconcile against the earlier cohort")
    if not prior_path or not os.path.exists(prior_path):
        LOG.say("5.0 no prior cohort file supplied; reconciliation NOT RUN")
        LOG.decision("5.0 reconciliation", "not run; --prior not supplied")
        return None

    prior = []
    with open(prior_path, "r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            prior.append(row)
    prior_bs = {r["biosample"].strip() for r in prior if r.get("biosample")}
    prior_by_bs = {r["biosample"].strip(): r for r in prior if r.get("biosample")}
    LOG.count("5.1 isolates in the earlier cohort", len(prior_bs))

    new_by_bs = {r["biosample"]: r for r in african if r.get("biosample")}
    new_bs = set(new_by_bs)
    LOG.count("5.1 isolates in the new cohort", len(new_bs))

    retained = sorted(prior_bs & new_bs)
    added = sorted(new_bs - prior_bs)
    dropped = sorted(prior_bs - new_bs)
    LOG.count("5.1 retained", len(retained))
    LOG.count("5.1 newly added", len(added))
    LOG.count("5.1 dropped", len(dropped))

    cutoff = prior_query_date
    add_rows = []
    for bs in added:
        r = new_by_bs[bs]
        sub = (r.get("submission_date") or "").strip()[:10]
        if sub and cutoff and sub > cutoff:
            reason = "deposited after the earlier query date (%s > %s)" % (sub, cutoff)
        elif r.get("location_fallback_flag", "no") != "no":
            reason = "location absent from the Pathogen Detection location field; assigned by %s fallback" \
                     % r["location_fallback_flag"].replace("yes-", "")
        else:
            reason = "location string not matched by the earlier 45-name wildcard list"
        add_rows.append({
            "biosample": bs, "assigned_country": r.get("assigned_country", ""),
            "location_string": r.get("location_string", ""),
            "submission_date": sub, "mcr_variants": r.get("mcr_variants", ""),
            "reason_earlier_search_missed_it": reason,
        })

    drop_rows = []
    for bs in dropped:
        p = prior_by_bs[bs]
        drop_rows.append({
            "biosample": bs, "prior_country": p.get("country", ""),
            "prior_mcr": p.get("mcr", ""), "prior_pd_target": p.get("pd_target", ""),
            "reason_dropped": ("not present in build %s under the current selection rules; "
                               "check whether NCBI revised, re-called or suppressed the record"
                               % args.build_for_report),
        })

    found = {c: 0 for c in KNOWN_MISSING_COUNTRIES}
    for r in african:
        if r["assigned_country"] in found and r["biosample"] not in prior_bs:
            found[r["assigned_country"]] += 1
    for c, n in found.items():
        LOG.count("5.2 newly added isolates from %s" % c, n)
    total_known = sum(found.values())
    LOG.count("5.2 newly added isolates from the three known-missing countries", total_known)
    if total_known < 9:
        LOG.say("5.2 STOP CONDITION: fewer than nine isolates recovered from South Africa, the Central "
                "African Republic and the Democratic Republic of the Congo (found %d). Diagnose before "
                "continuing: check the residual list, the qualifier ruling in 2.3 and the alias table." % total_known)
        LOG.decision("5.2 stop condition", "TRIGGERED (%d of the expected 9 recovered)" % total_known)
    else:
        LOG.decision("5.2 stop condition", "not triggered (%d recovered from the three countries)" % total_known)

    return {"retained": retained, "added": add_rows, "dropped": drop_rows, "known_missing": found}


# --------------------------------------------------------------------------
# STEP 6 - outputs
# --------------------------------------------------------------------------

FINAL_FIELDS = [
    "pd_target", "biosample", "assembly", "sra_runs", "sra_platform", "sra_layout",
    "sra_total_spots", "sra_total_bases", "sra_source", "organism",
    "submitted_species_is_ecoli", "species_decision", "strain", "host", "isolation_source",
    "collection_date", "location_string", "bioproject", "submission_date", "epi_type",
    "amr_genotype", "mcr_variants", "mcr_qualifiers", "mcr_has_partial",
    "mcr_has_internal_stop", "mcr_excluded_calls", "assigned_country", "country_match_alias",
    "country_assignment_reason", "location_fallback_flag", "western_sahara_review",
    "assigned_sector", "lat_lon",
]

GLOBAL_FIELDS = [
    "pd_target", "biosample", "assembly", "organism", "taxgroup", "strain", "host",
    "isolation_source", "collection_date", "location_string", "lat_lon", "bioproject",
    "submission_date", "epi_type", "amr_genotype", "mcr_variants", "mcr_qualifiers",
    "mcr_has_partial", "mcr_has_internal_stop", "mcr_excluded_calls",
]


def write_tsv(path, fields, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (r.get(k, "") if r.get(k) is not None else "") for k in fields})
    LOG.say("wrote %s (%d rows)" % (path, len(rows)))


def main():
    ap = argparse.ArgumentParser(description="Rebuild the mcr-positive African E. coli / Shigella cohort")
    ap.add_argument("--outdir", default="mcr_rebuild")
    ap.add_argument("--build", default=None, help="pin a PDG build, e.g. PDG000000004.6345")
    ap.add_argument("--prior", default=None, help="TSV of the earlier cohort with a 'biosample' column")
    ap.add_argument("--prior-date", default="2026-09-26", help="date the earlier search was run (YYYY-MM-DD)")
    ap.add_argument("--email", default=None, help="contact address sent to NCBI, strongly recommended")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    ap.add_argument("--pause", type=float, default=None, help="seconds between E-utilities calls")
    ap.add_argument("--skip-eutils", action="store_true")
    ap.add_argument("--offline", action="store_true",
                    help="no network at all; requires --build and --reuse-metadata, implies --skip-eutils")
    ap.add_argument("--reuse-metadata", default=None, help="path to an already-downloaded metadata TSV")
    ap.add_argument("--keep", nargs="*", default=None, help="qualifiers to keep in 2.3")
    ap.add_argument("--drop", nargs="*", default=None, help="qualifiers to exclude in 2.3")
    args = ap.parse_args()

    if args.offline:
        args.skip_eutils = True
    if args.pause is None:
        args.pause = 0.11 if args.api_key else 0.35
    os.makedirs(args.outdir, exist_ok=True)

    if not args.email:
        LOG.say("WARNING no --email supplied. NCBI asks for a contact address with E-utilities traffic.")

    meta = step1(args)
    args.build_for_report = meta["build"]

    rows, cols, header, versions = step2(meta, args)
    african, residual, alias_hits, excl_hits, tests = step3(rows, args)
    collisions, to_screen, query_date = step4(african, args)
    recon = step5(african, args.prior, args.prior_date, args)

    od = args.outdir
    # 6.1
    write_tsv(os.path.join(od, "01_final_cohort.tsv"), FINAL_FIELDS,
              sorted(african, key=lambda r: (r["assigned_country"], r["biosample"])))
    # 6.2
    write_tsv(os.path.join(od, "02_global_mcr_positive.tsv"), GLOBAL_FIELDS, rows)
    # 6.3
    write_tsv(os.path.join(od, "03_unresolved_location_residual.tsv"),
              ["pd_target", "biosample", "assembly", "organism", "location_string", "lat_lon",
               "host", "isolation_source", "collection_date", "bioproject", "submission_date",
               "mcr_variants", "country_assignment_reason"], residual)
    # 6.4 alias table with per-term match counts, including terms that matched none
    alias_rows = []
    for canon, (a2, a3, aliases) in AU_STATES.items():
        emitted = set()
        for al in aliases:
            if al in emitted:
                continue
            emitted.add(al)
            alias_rows.append({"canonical_country": canon, "iso_alpha2": a2, "iso_alpha3": a3,
                               "alias_term": al, "normalised_term": norm(al), "term_type": "alias",
                               "isolates_matched": alias_hits.get((canon, al), 0)})
        for code, ttype in ((a3, "iso_alpha3"), (a2, "iso_alpha2")):
            if code in emitted:
                continue
            emitted.add(code)
            alias_rows.append({"canonical_country": canon, "iso_alpha2": a2, "iso_alpha3": a3,
                               "alias_term": code, "normalised_term": norm(code), "term_type": ttype,
                               "isolates_matched": alias_hits.get((canon, code), 0)})
    for terr, names in NON_AU_EXCLUSIONS.items():
        for nme in names:
            alias_rows.append({"canonical_country": "(excluded non-AU territory) " + terr,
                               "iso_alpha2": "", "iso_alpha3": "", "alias_term": nme,
                               "normalised_term": norm(nme), "term_type": "exclusion",
                               "isolates_matched": excl_hits.get(terr, 0) if nme == names[0] else 0})
    write_tsv(os.path.join(od, "04_country_alias_lookup.tsv"),
              ["canonical_country", "iso_alpha2", "iso_alpha3", "alias_term", "normalised_term",
               "term_type", "isolates_matched"], alias_rows)
    LOG.count("6.4 alias terms in the lookup table", len(alias_rows))
    LOG.count("6.4 alias terms that matched nothing", sum(1 for a in alias_rows if not a["isolates_matched"]))

    # species screen list
    write_tsv(os.path.join(od, "05_species_screen_list.tsv"),
              ["biosample", "assembly", "organism", "assigned_country", "mcr_variants", "species_decision"],
              to_screen)

    # reconciliation
    if recon:
        write_tsv(os.path.join(od, "06_reconciliation_added.tsv"),
                  ["biosample", "assigned_country", "location_string", "submission_date",
                   "mcr_variants", "reason_earlier_search_missed_it"], recon["added"])
        write_tsv(os.path.join(od, "07_reconciliation_dropped.tsv"),
                  ["biosample", "prior_country", "prior_mcr", "prior_pd_target", "reason_dropped"],
                  recon["dropped"])
        with open(os.path.join(od, "08_reconciliation_retained.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(recon["retained"]) + "\n")

    # 6.6 final counts
    by_country = Counter(r["assigned_country"] for r in african)
    by_variant = Counter()
    for r in african:
        for v in (r["mcr_variants"] or "").split(";"):
            if v:
                by_variant[v] += 1
    with open(os.path.join(od, "09_final_counts.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "key", "value"])
        w.writerow(["isolates_retained", "", len(african)])
        w.writerow(["countries_represented", "", len(by_country)])
        for k, v in sorted(by_country.items()):
            w.writerow(["by_country", k, v])
        for k, v in sorted(by_variant.items()):
            w.writerow(["by_mcr_variant", k, v])
        for k, v in sorted(Counter(r["assigned_sector"] for r in african).items()):
            w.writerow(["by_sector", k, v])
    LOG.count("6.6 isolates retained", len(african))
    LOG.count("6.6 countries represented", len(by_country))

    # boundary test report + collisions into the log
    tb = "\n".join("%-6s %-46s expected=%-34s got=%s" % (st, loc, exp, got)
                   for (loc, exp, got, _al, _rs, st) in tests)
    cb = "\n".join("%-18s %-24s %s" % c for c in collisions) or "(none)"
    vb = "\n".join("%-28s %s" % (k, "; ".join("%s (n=%d)" % (x, n) for x, n in v)) for k, v in versions.items()) \
         or "(no version column found in the metadata table)"
    prov = ("build accession      %s\nbuild directory      %s\nmetadata file        %s\n"
            "metadata url         %s\nmetadata bytes       %d\nmetadata sha256      %s\n"
            "metadata md5         %s\nLast-Modified        %s\nE-utilities date     %s\n"
            "earlier query date   %s\nscript version       %s\npython               %s\n"
            % (meta["build"], meta["build_url"], meta["metadata_file"], meta["metadata_url"],
               meta["metadata_bytes"], meta["metadata_sha256"], meta["metadata_md5"],
               meta["metadata_last_modified"], query_date, args.prior_date, VERSION,
               sys.version.replace("\n", " ")))

    LOG.dump(os.path.join(od, "10_log.txt"), extra_sections=[
        ("PROVENANCE", prov),
        ("STEP 1.4 VERSIONS AS STATED IN THE METADATA TABLE", vb),
        ("STEP 2.3 QUALIFIER RULING AND RATIONALE", STEP23_RATIONALE),
        ("STEP 3.4 BOUNDARY TEST REPORT", tb),
        ("STEP 4.3 DEDUPLICATION COLLISIONS", cb),
        ("BUILD DIRECTORY LISTING", ", ".join(meta["build_dir_entries"])),
        ("METADATA DIRECTORY LISTING", ", ".join(meta["metadata_dir_entries"])),
        ("METADATA TABLE HEADER", ", ".join(header)),
    ])
    LOG.say("done. outputs in %s" % os.path.abspath(od))


if __name__ == "__main__":
    main()
