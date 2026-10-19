# Repository Guidelines

## Project Structure

This is the Donwulff/Exomiser fork, a Java 21 Maven project. Read `README.md`
for fork behaviour; `docs/` contains reStructuredText user documentation.

- `exomiser-core`: analysis, filters, prioritisers, variant models, DAOs, and writers.
- `exomiser-spring-boot-autoconfigure`: data paths and Spring wiring; starter/test
  modules support embedding.
- `exomiser-cli`, `exomiser-web`, `exomiser-rest-prioritiser`: application entry points.
- `exomiser-data-genome` and `exomiser-data-phenotype`: build reference bundles.
- Modules use `src/main/java`, `src/test/java`, and `src/*/resources`.
  CLI examples are in `exomiser-cli/src/main/resources/examples/`.

Pathogenicity changes span core, autoconfiguration, and writers.
Core HTML templates are `results_bootstrap_5.html` (current) and `results.html`
(legacy); web has its own `results.html`.
Variant links use core `report-links.yml`.

## Build, Test, and Run

Use `./mvnw` with JDK 21. Check `./mvnw -version`; set `JAVA_HOME` if necessary.
Newer JDKs can fail Mockito/ByteBuddy tests.

```bash
./mvnw clean install -DskipTests=true          # Build/install all modules
./mvnw -pl exomiser-cli -am package -DskipTests=true  # Package CLI and dependencies
./mvnw test -pl exomiser-core                  # Core tests
./mvnw test -pl exomiser-core -Dtest=ReportLinksTest  # Focused regression
./mvnw test                                   # Full suite
java -Xmx4g -jar exomiser-cli/target/exomiser-cli-15.1.0.jar analyse \
  --analysis /path/to/analysis.yml --spring.config.location=file:./application.properties
```

Tests use JUnit Jupiter, Mockito, and Hamcrest. Name classes `*Test`; add
regressions using synthetic/public fixtures. Never enable parallel Surefire:
MVStore requires exclusive access. JaCoCo has no enforced coverage minimum.

## Coding Conventions

Follow surrounding style: four-space indentation, `UpperCamelCase` classes,
`lowerCamelCase` members, and `UPPER_SNAKE_CASE` constants. Preserve AGPL headers.
No global formatter/checker is configured.

## Data and Configuration

Never commit real VCF/BAM/CRAM/PED files, patient configurations, reports, logs,
reference ZIPs, or databases. `application.properties`, `cases/`, `results/`,
and `logs/` at the root are gitignored. The properties template is in CLI resources.
Ignore misplaced payloads; never stage them.

Data lives separately, locally under `/mnt/GenomicData/Exomiser`; consult its
`AGENTS.md` before operating there. Never search that volume recursively without
bounds. `~/Exomiser` is an old installation, not this build tree.

`data/` holds small torrents and canonical `seed-torrents.sh` only. Torrent official
release ZIPs with size-matching BEP-19 webseeds, never extracted trees:

```bash
./data/seed-torrents.sh check
./data/seed-torrents.sh create --release 2602
```

`EXOMISER_DATA_DIR` selects the payload directory; edit `EXTRA_GLOBS` for optional
payloads, kept outside Git. Data releases can serve multiple software releases.

## Commits and Pull Requests

Use imperative subjects, e.g. “Add configurable report links”. Describe
behaviour, relevant issues, and validation in PRs. Use synthetic data in screenshots.
Review staged files for sensitive data. Preserve unrelated working changes.

This is the shared guide; `CLAUDE.md` imports it.
`IDEAS.md` and session exports are local historical context, not an approved
roadmap or files to commit automatically.
