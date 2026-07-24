#!/usr/bin/env bash
#
# seed-torrents.sh — create and/or seed BitTorrent files for the public Exomiser
# data bundles in this directory.
#
# CANONICAL HOME in git: data/seed-torrents.sh in the Donwulff/Exomiser fork
# (https://github.com/Donwulff/Exomiser), next to published *.torrent files.
# When this file lives on the data volume (/mnt/GenomicData/Exomiser or similar)
# it is the LOCAL OPERATIONAL copy — run it where the zip payloads and seeder
# live. Keep both in sync; edit the repo copy as the source of truth.
# See repo AGENTS.md and the data-volume AGENTS.md.
#
# Background: these .zip bundles are public data. We seed them over BitTorrent
# because direct HTTP downloads were being corrupted in transit. As of 2026-07-24
# upstream distributes the bundles over an HTTP/Globus CDN only — no official
# .torrent files could be found — so:
#
#   1. MINT YOUR OWN (the working path today) — `create` runs mktorrent over the
#      zips with public trackers + an HTTP webseed pointing at the upstream CDN
#      (WEBSEED_BASE below). The webseed lets downloaders fall back to the CDN when
#      no peer is seeding, so the torrent works even with a cold swarm.
#
#   2. SEED UPSTREAM'S TORRENTS — only if/when upstream ever publishes .torrent
#      files: add the official .torrent to a client pointed at THIS directory and
#      `--verify`; a matching infohash seeds the existing swarm. See `seed_help`.
#
# Torrent the .zip files, NOT the extracted directories: the zip is the
# distributable artifact and is what any upstream infohash refers to.
#
# Usage:
#   ./seed-torrents.sh check           # report tooling + which files are present
#   ./seed-torrents.sh create          # mktorrent for each FILES entry -> *.torrent
#   ./seed-torrents.sh show            # print infohash + magnet for existing *.torrent
#   ./seed-torrents.sh seed_help       # how to seed (upstream or self-made) via transmission
#
# Requires (install yourself; not auto-installed):
#   sudo apt install mktorrent transmission-daemon transmission-cli
#
set -euo pipefail
cd "$(dirname "$0")"

# ---- config -----------------------------------------------------------------

# Bundles to seed. Add the NCBoost file here if you also want to seed it.
FILES=(
  "2602_hg19.zip"
  "2602_hg38.zip"
  "2602_phenotype.zip"
  # "ncboost_score_hg19_v20190108.tsv.gz"
)

# Public trackers (UDP). Trim/extend as you like.
TRACKERS=(
  "udp://tracker.opentrackr.org:1337/announce"
  "udp://open.stealth.si:80/announce"
  "udp://tracker.torrent.eu.org:451/announce"
  "udp://exodus.desync.com:6969/announce"
  "udp://tracker.openbittorrent.com:6969/announce"
)

# HTTP(S) webseed base (BEP-19). If set, downloaders can pull from the original
# server when no peer is seeding — hugely helpful for a cold swarm. The full URL
# used per file is  ${WEBSEED_BASE}<filename> , so it MUST directly serve e.g.
#   ${WEBSEED_BASE}2602_hg19.zip
# Below is the official 2602 host (Globus CDN from the Exomiser install docs),
# verified 2026-07-24 to return HTTP 200 with Content-Length byte-matching the
# local zips. NOTE: Globus guest-collection hostnames can rotate between releases
# — re-verify with `curl -sI "${WEBSEED_BASE}2602_hg19.zip"` before regenerating.
WEBSEED_BASE="https://g-879a9f.f5dc97.75bc.dn.glob.us/data/"

# Piece size exponent: piece = 2^N bytes. 24 = 16 MiB, sensible for ~20 GB files
# (keeps the .torrent small). Omit -l to let mktorrent auto-pick.
PIECE_EXP=24

# ---- helpers ----------------------------------------------------------------

have() { command -v "$1" >/dev/null 2>&1; }

check() {
  echo "== tooling =="
  for t in mktorrent transmission-show transmission-remote transmission-daemon; do
    if have "$t"; then echo "  ok   $t ($(command -v "$t"))"; else echo "  MISS $t"; fi
  done
  echo "== files =="
  for f in "${FILES[@]}"; do
    if [[ -f "$f" ]]; then
      printf "  ok   %-40s %s\n" "$f" "$(du -h "$f" | cut -f1)"
    else
      printf "  MISS %-40s (not found)\n" "$f"
    fi
  done
  if [[ -z "$WEBSEED_BASE" ]]; then
    echo "note: WEBSEED_BASE is empty — torrents will have no HTTP webseed."
  else
    echo "webseed: $WEBSEED_BASE"
  fi
}

create() {
  have mktorrent || { echo "mktorrent not installed: sudo apt install mktorrent" >&2; exit 1; }
  local announce=() ; for t in "${TRACKERS[@]}"; do announce+=(-a "$t"); done
  for f in "${FILES[@]}"; do
    [[ -f "$f" ]] || { echo "skip (missing): $f" >&2; continue; }
    local out="${f}.torrent"
    if [[ -e "$out" ]]; then echo "skip (exists): $out — remove it to regenerate" >&2; continue; fi
    local web=()
    [[ -n "$WEBSEED_BASE" ]] && web=(-w "${WEBSEED_BASE}${f}")
    echo ">> mktorrent $f -> $out"
    mktorrent -l "$PIECE_EXP" "${announce[@]}" "${web[@]}" -o "$out" "$f"
  done
  echo "done. Verify with: $0 show"
}

show() {
  have transmission-show || { echo "transmission-show not installed: sudo apt install transmission-cli" >&2; exit 1; }
  shopt -s nullglob
  local any=0
  for tf in *.torrent; do
    any=1
    echo "== $tf =="
    transmission-show "$tf" | grep -E "Name:|Hash:|Total Size:|Piece Count:|Magnet:" || transmission-show "$tf"
    echo
  done
  [[ $any -eq 0 ]] && echo "no *.torrent files here yet — run: $0 create"
}

seed_help() {
  cat <<'EOF'
Seeding with transmission-daemon (headless), pointed at THIS directory so the
existing zips are used as-is:

  # 1. install + start the daemon (once)
  sudo apt install transmission-daemon transmission-cli
  sudo systemctl stop transmission-daemon      # edit config while stopped
  #    set "download-dir" to this directory in
  #    /etc/transmission-daemon/settings.json  (or run as your user with --download-dir)
  sudo systemctl start transmission-daemon

  # 2a. SEED UPSTREAM'S TORRENT — download their .torrent for each zip, then:
  transmission-remote -a 2602_hg19.zip.torrent --download-dir "$PWD"
  transmission-remote -t all --verify          # rechecks local data -> seeds if it matches

  # 2b. SEED YOUR OWN — after `./seed-torrents.sh create`:
  transmission-remote -a 2602_hg19.zip.torrent --download-dir "$PWD"
  transmission-remote -t all --verify

  # check status
  transmission-remote -l

Key point: --verify makes transmission hash the zips already on disk; if they
match the torrent's pieces it skips downloading and starts seeding immediately.
Do not rename/recompress a zip you are seeding — that breaks the infohash match.
EOF
}

case "${1:-check}" in
  check)     check ;;
  create)    create ;;
  show)      show ;;
  seed_help) seed_help ;;
  *) echo "usage: $0 {check|create|show|seed_help}" >&2; exit 2 ;;
esac
