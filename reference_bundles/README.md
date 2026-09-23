# Adding a plant reference

Plant- and assembly-specific data belongs in a reference bundle. Shared files
under `scripts/` should not contain chromosome names, chromosome lengths,
organism database packages, or coordinates.

## UI reference wizard

Most users do not need to write YAML. `Methylome.Plants_UI.sh` collects the organism, assembly, GTF/GFF3 and optional resources, validates them, and generates a reusable bundle under `reference_bundles/generated/`.

GO/OrgDb and KEGG are selected independently. When the chosen OrgDb provides a `CHRLENGTHS` object, the wizard can use it as a chromosome-length source after checking its sequence names and coordinates against the annotation.

The manual workflow below remains useful for version-controlled reference packages.

## 1. Create the bundle

Copy `template.yaml` to a species/assembly-named YAML file. Keep supporting
files in a sibling directory so every relative path remains portable.

Required fields:

- `schema_version`
- `species.id`, `species.display_name`, and `species.assembly`
- `genome.primary_seqlevels`
- `annotation.genes`

Exact chromosome lengths are strongly recommended but are not mandatory for
core analysis. The preferred sources, in order, are a reference FASTA/FAI, a
two-column chromosome-size table, a compatible BSgenome, or `CHRLENGTHS` from
the selected OrgDb. When supplied as a table, use `seqname` and `length`
columns. The names in `primary_seqlevels` are the canonical names used in
pipeline output.

## 2. Declare sequence aliases

Add a tab-delimited `seqname_aliases` file with `alias` and `canonical`
columns. List only explicit, unambiguous mappings, for example:

```text
alias	canonical
1	Chr01
chr01	Chr01
CP	Chloroplast
```

Do not encode aliases with regular-expression guesses. Two input sequence
levels may not collapse to the same canonical name within one input object.

## 3. Configure annotations

Gene annotations may be GFF3, GTF, or CSV. Map source metadata fields with the
`annotation.*_field` settings. GFF3 `ID`/`Parent` ancestry is resolved
automatically and the source attributes are preserved.

TE annotations may be BED, GFF3/GTF, CSV, or TSV. Coordinates and a TE ID are
required. Family and superfamily are optional; absent classifications are
reported as unclassified.

Descriptions, organelle names, centromeres, heterochromatin, TFBS, GO, KEGG,
gbM sets, and functional gene sets are optional. Analyses that need an absent
resource are disabled with a message; core methylation and DMR analyses remain
available.

## 4. Keep assemblies consistent

The methylation calls, FASTA used for alignment, gene annotation, TE
annotation, coordinate tracks, and chromosome-size table must describe the
same assembly. Bundle BED files use zero-based, half-open coordinates.

At startup the pipeline:

1. resolves bundle paths relative to the YAML file;
2. validates configured files and chromosome sizes;
3. applies only declared sequence aliases;
4. filters to configured primary sequences when requested;
5. checks sequence overlap and coordinate bounds; and
6. writes `reference_bundle_resolved.yaml` into the comparison results.

The UI performs an earlier validation pass and also checks GTF gene IDs against
the selected OrgDb key types. A KEGG organism code alone does not guarantee
that GTF identifiers match KEGG gene identifiers; add a `gene_id`/`kegg_id`
mapping table when needed.

## 5. Validate before a full run

Run the repository smoke test:

```bash
Rscript scripts/test_scripts/reference_bundle_smoke.R
Rscript scripts/test_scripts/reference_wizard_smoke.R
```

Then run a small comparison with the new bundle:

```bash
./scripts/Methylome.Plants.sh samples.txt \
  --reference_bundle reference_bundles/genus_species_assembly.yaml
```

Use `--annotation_file`, `--description_file`, or `--TEs_file` only for
temporary overrides. Permanent species support should stay in the bundle.
