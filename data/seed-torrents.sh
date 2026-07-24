#!/usr/bin/env bash
#
# seed-torrents.sh — mint BitTorrent metainfo for public Exomiser data zips.
#
# Canonical home in git: data/seed-torrents.sh (Donwulff/Exomiser).
# Zips live on a separate data volume; only tiny *.torrent files belong in git.
#
# Automagic flow (default when you have zips on disk):
#   1. Discover YYMM_{hg19,hg38,phenotype}.zip under DATA_DIR
#   2. Probe known HTTP bases until Content-Length matches local size
#      (BEP-19 webseed = official zip fallback for a cold swarm)
#   3. mktorrent → TORRENT_DIR (this directory in the repo)
#
# Never torrents extracted trees — only the official .zip artifacts.
#
# Usage:
#   ./seed-torrents.sh                     # same as check
#   ./seed-torrents.sh check               # tooling + discovered zips + webseed probe
#   ./seed-torrents.sh create              # mint missing *.torrent into TORRENT_DIR
#   ./seed-torrents.sh create --force      # regenerate even if .torrent exists
#   ./seed-torrents.sh create --release 2602
#   ./seed-torrents.sh verify-webseed      # HEAD-check webseed sizes only
#   ./seed-torrents.sh show                # infohash / magnet for repo torrents
#   ./seed-torrents.sh seed_help           # transmission-daemon notes
#
# Environment / flags (flags win):
#   EXOMISER_DATA_DIR / --data-dir DIR     zip location (default: /mnt/GenomicData/Exomiser)
#   EXOMISER_TORRENT_DIR / --torrent-dir   where to write *.torrent (default: this script's dir)
#   EXOMISER_WEBSEED_BASE / --webseed URL  force webseed base (must end with /)
#   EXOMISER_RELEASE / --release YYMM      only that release (else all complete sets + extras)
#
# Requires: mktorrent; curl for webseed probe; transmission-cli optional for show.
#   sudo apt install mktorrent curl transmission-cli
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Prefer system mktorrent; fall back to a dpkg-extracted copy under data/.tools/
if ! command -v mktorrent >/dev/null 2>&1; then
  _local_mk="$SCRIPT_DIR/.tools/mktorrent-pkg/usr/bin"
  [[ -x "$_local_mk/mktorrent" ]] && PATH="$_local_mk:$PATH"
fi
unset _local_mk

# ---- defaults ---------------------------------------------------------------

DATA_DIR="${EXOMISER_DATA_DIR:-/mnt/GenomicData/Exomiser}"
TORRENT_DIR="${EXOMISER_TORRENT_DIR:-$SCRIPT_DIR}"
WEBSEED_BASE="${EXOMISER_WEBSEED_BASE:-}"
RELEASE="${EXOMISER_RELEASE:-}"
FORCE=0
# Optional extras (basename patterns under DATA_DIR) — not part of a YYMM trio
EXTRA_GLOBS=(
  # "ncboost_score_*.tsv.gz"
)

# Known official HTTP bases (tried in order until sizes match). Globus guest
# collection hostnames can rotate per release — probe, don't hard-require one.
WEBSEED_CANDIDATES=(
  "https://g-879a9f.f5dc97.75bc.dn.glob.us/data/"
  "https://data.monarchinitiative.org/exomiser/latest/"
)

TRACKERS=(
  "udp://tracker.opentrackr.org:1337/announce"
  "udp://open.stealth.si:80/announce"
  "udp://tracker.torrent.eu.org:451/announce"
  "udp://exodus.desync.com:6969/announce"
  "udp://tracker.openbittorrent.com:6969/announce"
)

# piece = 2^N bytes; 24 = 16 MiB — keeps .torrent small for ~20 GB files
PIECE_EXP=24

BUNDLE_SUFFIXES=(hg19 hg38 phenotype)

# ---- arg parse --------------------------------------------------------------

CMD="${1:-check}"
if [[ $# -gt 0 ]]; then shift; fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir)      DATA_DIR="$2"; shift 2 ;;
    --torrent-dir)   TORRENT_DIR="$2"; shift 2 ;;
    --webseed)       WEBSEED_BASE="$2"; shift 2 ;;
    --release)       RELEASE="$2"; shift 2 ;;
    --force|-f)      FORCE=1; shift ;;
    -h|--help)       CMD=help; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

# normalize trailing slash on webseed if user set one
if [[ -n "$WEBSEED_BASE" && "$WEBSEED_BASE" != */ ]]; then
  WEBSEED_BASE="${WEBSEED_BASE}/"
fi

# ---- helpers ----------------------------------------------------------------

have() { command -v "$1" >/dev/null 2>&1; }

die() { echo "error: $*" >&2; exit 1; }

usage() {
  sed -n '2,35p' "$0" | sed 's/^# \?//'
}

# Complete release = all three zips present. Echo YYMM tags, sorted.
discover_releases() {
  local d="$DATA_DIR" s tag
  [[ -d "$d" ]] || return 0
  shopt -s nullglob
  local -A seen=()
  for s in "${BUNDLE_SUFFIXES[@]}"; do
    for f in "$d"/[0-9][0-9][0-9][0-9]_"$s".zip; do
      tag="$(basename "$f")"
      tag="${tag%%_*}"
      seen["$tag"]=1
    done
  done
  local t complete
  for t in $(printf '%s\n' "${!seen[@]}" | sort); do
    complete=1
    for s in "${BUNDLE_SUFFIXES[@]}"; do
      [[ -f "$d/${t}_${s}.zip" ]] || complete=0
    done
    [[ $complete -eq 1 ]] && echo "$t"
  done
}

# Print basenames of zips we should torrent for the selected release(s).
list_target_zips() {
  local d="$DATA_DIR" t s f
  if [[ -n "$RELEASE" ]]; then
    for s in "${BUNDLE_SUFFIXES[@]}"; do
      f="${RELEASE}_${s}.zip"
      [[ -f "$d/$f" ]] || die "missing $d/$f (release $RELEASE)"
      echo "$f"
    done
  else
    local any=0
    while IFS= read -r t; do
      any=1
      for s in "${BUNDLE_SUFFIXES[@]}"; do
        echo "${t}_${s}.zip"
      done
    done < <(discover_releases)
    # partial latest? if no complete set, still list any matching YYMM zips
    if [[ $any -eq 0 ]]; then
      shopt -s nullglob
      for f in "$d"/[0-9][0-9][0-9][0-9]_hg19.zip \
               "$d"/[0-9][0-9][0-9][0-9]_hg38.zip \
               "$d"/[0-9][0-9][0-9][0-9]_phenotype.zip; do
        basename "$f"
      done
    fi
  fi
  # extras
  local g
  shopt -s nullglob
  for g in "${EXTRA_GLOBS[@]}"; do
    for f in "$d"/$g; do
      [[ -f "$f" ]] && basename "$f"
    done
  done
}

local_size() {
  stat -c '%s' "$1" 2>/dev/null || stat -f '%z' "$1"
}

# HEAD Content-Length for URL; empty if missing/failed.
remote_length() {
  local url="$1" len
  len="$(curl -sI -L --max-time 30 "$url" 2>/dev/null \
    | tr -d '\r' \
    | awk 'tolower($1)=="content-length:" {print $2; exit}')"
  [[ -n "$len" && "$len" =~ ^[0-9]+$ ]] && echo "$len"
}

# Does base serve this basename with matching size? 0=yes.
webseed_matches() {
  local base="$1" name="$2" path="$3"
  local want got
  want="$(local_size "$path")"
  got="$(remote_length "${base}${name}")"
  [[ -n "$got" && "$got" == "$want" ]]
}

# Probe WEBSEED_BASE or candidates against first target zip; set WEBSEED_BASE.
resolve_webseed() {
  local first path name base
  first="$(list_target_zips | head -n1 || true)"
  [[ -n "$first" ]] || { echo "  (no local zips to probe against)"; return 1; }
  path="$DATA_DIR/$first"
  name="$first"

  if [[ -n "$WEBSEED_BASE" ]]; then
    if webseed_matches "$WEBSEED_BASE" "$name" "$path"; then
      echo "  ok   forced webseed matches $name ($(local_size "$path") bytes)"
      return 0
    fi
    echo "  FAIL forced webseed size mismatch or unreachable: ${WEBSEED_BASE}${name}" >&2
    echo "       local=$(local_size "$path") remote=$(remote_length "${WEBSEED_BASE}${name}" || echo none)" >&2
    return 1
  fi

  for base in "${WEBSEED_CANDIDATES[@]}"; do
    echo -n "  probe $base ... "
    if webseed_matches "$base" "$name" "$path"; then
      WEBSEED_BASE="$base"
      echo "MATCH $name"
      return 0
    fi
    echo "no"
  done
  echo "  WARN no webseed base matched local sizes — torrents will have no HTTP fallback" >&2
  WEBSEED_BASE=""
  return 1
}

# Range-request smoke test. Some CDNs ignore Range on HEAD; use a 1-byte GET.
range_ok() {
  local url="$1" hdr
  hdr="$(curl -sL --max-time 30 -r 0-0 -D - -o /dev/null "$url" 2>/dev/null | tr -d '\r')"
  grep -qiE 'HTTP/.* 206|Content-Range:' <<<"$hdr"
}

check() {
  echo "== tooling =="
  for t in mktorrent curl transmission-show transmission-remote; do
    if have "$t"; then echo "  ok   $t ($(command -v "$t"))"; else echo "  MISS $t"; fi
  done
  echo "== paths =="
  echo "  DATA_DIR    = $DATA_DIR"
  echo "  TORRENT_DIR = $TORRENT_DIR"
  [[ -d "$DATA_DIR" ]] || echo "  WARN DATA_DIR missing"
  mkdir -p "$TORRENT_DIR"

  echo "== releases on disk =="
  local r
  r="$(discover_releases | tr '\n' ' ')"
  if [[ -n "${r// /}" ]]; then
    echo "  complete: $r"
  else
    echo "  (no complete YYMM_hg19+hg38+phenotype set found)"
  fi
  [[ -n "$RELEASE" ]] && echo "  filter: --release $RELEASE"

  echo "== target zips =="
  local f any=0
  while IFS= read -r f; do
    any=1
    if [[ -f "$DATA_DIR/$f" ]]; then
      printf "  ok   %-40s %s\n" "$f" "$(du -h "$DATA_DIR/$f" | cut -f1)"
    else
      printf "  MISS %-40s\n" "$f"
    fi
  done < <(list_target_zips)
  [[ $any -eq 1 ]] || echo "  (none)"

  echo "== webseed =="
  if resolve_webseed; then
    local sample
    sample="$(list_target_zips | head -n1)"
    if range_ok "${WEBSEED_BASE}${sample}"; then
      echo "  ok   Accept-Ranges / 206 looks good for ${sample}"
    else
      echo "  WARN range request did not clearly return 206 for ${sample}"
    fi
    # verify every target
    while IFS= read -r f; do
      if webseed_matches "$WEBSEED_BASE" "$f" "$DATA_DIR/$f"; then
        echo "  ok   $f"
      else
        echo "  FAIL $f size mismatch vs ${WEBSEED_BASE}${f}" >&2
      fi
    done < <(list_target_zips)
  fi

  echo "== existing torrents in TORRENT_DIR =="
  shopt -s nullglob
  local tf
  for tf in "$TORRENT_DIR"/*.torrent; do
    echo "  $(basename "$tf")"
  done
  shopt -u nullglob
}

verify_webseed() {
  resolve_webseed || true
  [[ -n "$WEBSEED_BASE" ]] || die "no working webseed (set --webseed or fix CDN)"
  local f rc=0
  while IFS= read -r f; do
    if webseed_matches "$WEBSEED_BASE" "$f" "$DATA_DIR/$f"; then
      echo "ok  $f  <- ${WEBSEED_BASE}${f}"
    else
      echo "FAIL $f  local=$(local_size "$DATA_DIR/$f") remote=$(remote_length "${WEBSEED_BASE}${f}" || echo none)" >&2
      rc=1
    fi
  done < <(list_target_zips)
  return $rc
}

create() {
  have mktorrent || die "mktorrent not installed — sudo apt install mktorrent"
  have curl || die "curl required for webseed probe"
  [[ -d "$DATA_DIR" ]] || die "DATA_DIR not found: $DATA_DIR"
  mkdir -p "$TORRENT_DIR"

  local targets=()
  local f
  while IFS= read -r f; do targets+=("$f"); done < <(list_target_zips)
  [[ ${#targets[@]} -gt 0 ]] || die "no target zips under $DATA_DIR"

  # Best-effort webseed; allow create without it if probe fails
  resolve_webseed || true
  if [[ -n "$WEBSEED_BASE" ]]; then
    echo "webseed: $WEBSEED_BASE"
  else
    echo "webseed: (none)"
  fi

  local announce=() t
  for t in "${TRACKERS[@]}"; do announce+=(-a "$t"); done

  local out web path name
  for name in "${targets[@]}"; do
    path="$DATA_DIR/$name"
    [[ -f "$path" ]] || { echo "skip (missing): $path" >&2; continue; }
    out="$TORRENT_DIR/${name}.torrent"
    if [[ -e "$out" && $FORCE -eq 0 ]]; then
      echo "skip (exists): $out — pass --force to regenerate"
      continue
    fi
    [[ -e "$out" && $FORCE -eq 1 ]] && rm -f "$out"

    web=()
    [[ -n "$WEBSEED_BASE" ]] && web=(-w "${WEBSEED_BASE}${name}")

    # cd into DATA_DIR so the torrent "name" is the basename only
    echo ">> mktorrent $name -> $out  (hashing $(du -h "$path" | cut -f1); may take a while)"
    (
      cd "$DATA_DIR"
      mktorrent -l "$PIECE_EXP" "${announce[@]}" "${web[@]}" -o "$out" "$name"
    )
    echo "   wrote $out ($(du -h "$out" | cut -f1))"
  done
  echo "done. Review with: $0 show"
  echo "Commit only the new *.torrent files under data/ — never the zips."
}

show() {
  have transmission-show || die "transmission-show not installed — sudo apt install transmission-cli"
  shopt -s nullglob
  local tf any=0
  for tf in "$TORRENT_DIR"/*.torrent; do
    any=1
    echo "== $(basename "$tf") =="
    transmission-show "$tf" | grep -E "Name:|Hash:|Total Size:|Piece Count:|Magnet:|Webseed" \
      || transmission-show "$tf"
    echo
  done
  [[ $any -eq 1 ]] || echo "no *.torrent in $TORRENT_DIR — run: $0 create"
}

seed_help() {
  cat <<EOF
Seeding with transmission-daemon, pointed at DATA_DIR so existing zips are used:

  DATA_DIR=$DATA_DIR
  TORRENT_DIR=$TORRENT_DIR

  sudo apt install transmission-daemon transmission-cli mktorrent
  # set download-dir to DATA_DIR (settings.json or --download-dir)
  transmission-remote -a $TORRENT_DIR/2602_hg19.zip.torrent --download-dir "$DATA_DIR"
  transmission-remote -t all --verify
  transmission-remote -l

--verify hashes local zips against the torrent; matching pieces seed immediately.
Do not recompress/rename a zip you are seeding.

Regenerate torrents after a new data release lands in DATA_DIR:

  $0 create --release YYMM
  git add data/YYMM_*.zip.torrent && git commit -m "Torrents for Exomiser data YYMM"
EOF
}

case "$CMD" in
  check)          check ;;
  create)         create ;;
  verify-webseed) verify_webseed ;;
  show)           show ;;
  seed_help)      seed_help ;;
  help|-h|--help) usage ;;
  *) echo "usage: $0 {check|create|verify-webseed|show|seed_help} [options]" >&2; exit 2 ;;
esac
