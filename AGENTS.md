# Repository guidelines (agents & contributors)

This is the **Donwulff/Exomiser** software fork of
[exomiser/Exomiser](https://github.com/exomiser/Exomiser) — Java source, build,
tests, report templates, and small distribution helpers. It is **not** a place
to store genomes, VCFs, extracted data bundles, or analysis output.

## Keep real data and results out of this repository

**Do not commit** any of the following (and do not ask agents to add them):

| Kind | Examples | Why |
|------|----------|-----|
| Reference-data payloads | `*_hg19.zip`, extracted `YYMM_hg19/`, `*.mv.db`, large tabix score files | Multi‑GB public bundles; live on a dedicated data volume |
| Sequencing / clinical inputs | `*.vcf`, `*.vcf.gz`, BAM/CRAM, PED with real samples | Often sensitive; never belong in git |
| Case / analysis configs with real samples | Ad‑hoc YAML under `cases/`, patient HPO lists, local paths | Machine- and person-specific |
| Run output | HTML/TSV/JSON/VCF reports, `results/`, `logs/` | Regenerable; huge; may embed sample detail |
| Machine-local runtime config | Root `application.properties` with absolute data paths | Template is tracked under `exomiser-cli/src/main/resources/` |

What **is** fine in-repo:

- Source, tests, docs, example analysis YAMLs under
  `exomiser-cli/src/main/resources/examples/`
- **BitTorrent helpers only**: small `*.torrent` files and scripts under `data/`
  (see below) — not the zip payloads they describe
- Gitignored local paths for day-to-day runs (see root `.gitignore`)

`.gitignore` already excludes `/application.properties`, `/logs/`, `/results/`,
and `/cases/`. If something large or sample-specific still shows up in
`git status`, treat that as a bug: ignore it, do not stage it.

### Where data and results actually live

On this machine (illustrative; not part of the git tree):

- **Data volume / seeder host:** `/mnt/GenomicData/Exomiser` — release zips,
  extracted bundles, optional historical results, operational seeder.
  Agent notes for *that* directory: `/mnt/GenomicData/Exomiser/AGENTS.md`.
- **Runtime props:** repo-root `application.properties` (gitignored) →
  `exomiser.data-directory=…`.
- **Old install / case YAML reference only:** `~/Exomiser` — not this build tree.

Upstream (and this fork) expect a **separate** data directory; one data release
can serve many software point releases.

## Data distribution: torrents + official zip webseeds

The original fork shipped **`.torrent` files** for the public Exomiser data
bundles (HTTP downloads were unreliable/corrupt on some networks), with a
**BEP-19 webseed** back to the official zip URLs so a cold swarm still works.

Policy for continuing that:

1. **In this repo (`data/`):** version-controlled `*.torrent` files, and the
   canonical `data/seed-torrents.sh` used to mint them for new `YYMM` releases
   (update `FILES` / `WEBSEED_BASE`, re-verify Range support on the CDN host).
2. **On the data volume:** the actual `.zip` bytes and the running BitTorrent
   client. Keep an operational copy of the script next to the zips; edit the
   repo copy first, then mirror.
3. **Torrent the `.zip`s only — never the extracted trees.** There is no
   official fallback for unpacked directories, infohashes would not match any
   upstream artifact, and serving extracted data adds size without a real
   distribution benefit. Extract locally after download.
4. Optional extra files (e.g. NCBoost score tabix) may get torrents the same
   way; still do not commit the payload.

Do not commit extracted release directories or multi‑GB zips “for convenience”.

## Agent / working notes

- Prefer **this file** for durable “what belongs in the repo / what not to
  commit” rules.
- Deeper local architecture and build notes may live in untracked working
  files (`CLAUDE.md`, `IDEAS.md`, etc.). Whether those become tracked is a
  separate decision — do not assume they should be committed.
- Never run unbounded `find`/`grep` over `/mnt/GenomicData` (multi‑TB volume).

## Build and test (short)

Use the Maven wrapper and a Mockito-friendly JDK (prefer **21** for tests):

```bash
./mvnw clean install -DskipTests=true
./mvnw -pl exomiser-cli -am package -DskipTests=true
./mvnw test -pl exomiser-core
```

Do not enable parallel Surefire (MVStore single-open constraint). See module
docs and `README.md` for product behaviour; fork deltas (DANN/NCBOOST, report
links, multi-file CADD, call-quality UI) are documented at the top of
`README.md` / recent commit history.

## Commits

- Stage only source, tests, docs, and small helpers.
- Never stage `cases/`, `results/`, data zips/databases, or root runtime
  `application.properties`.
- Keep messages concise and focused on why the change exists.
