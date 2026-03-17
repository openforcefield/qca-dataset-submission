# GitHub Copilot Code Review Instructions

Use this document to review pull requests in this repository for **STANDARDS.md compliance** with high precision and low noise.

## Source of Truth and Precedence

1. `STANDARDS.md` is the primary authority.
2. This file explains how to apply `STANDARDS.md` in automated reviews.
3. If maintainers explicitly document an allowed exception in the PR discussion, do not keep reporting it as unresolved.

Start every review with: **“Experimental automated review; final authority is STANDARDS.md.”**

## Scope: What to Review

Review only structural/reproducibility/metadata compliance for dataset submissions.

Do review:
- Submission directory structure and naming under `submissions/`.
- Presence of required reproducibility files (generation script/notebook, environment file, README).
- Metadata completeness and consistency (name/version/descriptions/provenance/submitter/elements/charges/molecular weight stats).
- README completeness and consistency with dataset metadata.
- Versioning and revision/changelog rules.
- Root `README.md` table/index update in the correct dataset-type section.

Do not review:
- Scientific correctness of computed results.
- Notebook print verbosity, formatting style, or other non-standards nits.
- General code-quality opinions unless they directly block reproducibility or standards compliance.
- `dataset*.json.bz2` and `scaffold*.json.bz2` content-level review.

## Decision Rules (Critical to Avoid False Positives)

- Only mark **blocking** when a requirement is explicit in `STANDARDS.md` and currently unmet.
- If a requirement is expected “before merge” (for example DOI for force-field release datasets), phrase as: **“Must be present before merge; re-check at final review.”**
- Do not convert style suggestions into blockers.
- Do not repeat the same issue multiple times across files unless each location needs a distinct fix.
- If uncertain whether something is required, ask a single clarification question or mark as **manual verification needed**, not blocking.

## Required Checks

### 1) Submission Structure
- New dataset/revision in its own directory: `submissions/YYYY-MM-DD-{dataset-name-with-hyphens}`.
- Existing datasets are revised via version bump, not by mutating old released artifacts.
- Required files are present and described in local README.

### 2) Dataset Metadata (from `.py` / `.ipynb` setup code)
- Check metadata in generation script/notebook, not compressed data artifacts.
- Confirm required fields from `STANDARDS.md` are present.
- `long_description_url` should point to submission folder URL:
  `https://github.com/openforcefield/qca-dataset-submission/tree/master/submissions/{dir-name}`.
- Ensure dataset name uses spaces (not hyphens/underscores) in metadata.

### 3) README + Revision Policy
- Local submission README contains required metadata mirror and QC specification details.
- Changelog is required for revisions beyond initial submission; not required for initial submissions (e.g. `v4.0` or `v0.0`).
- README title should match dataset name.
- README file manifest should list/define files in the submission directory.

### 4) Versioning and Naming
- Version format must be `vX.Y`.
- Major version maps to standards conformance (`v4.x` for STANDARDS v4, `v0.x` for non-conforming/legacy cases).
- For `>=v4.0`, expected naming is `OpenFF {descriptive name} v{version}` with spaces.
- For force-field release datasets, expected naming uses `OpenFF SMIRNOFF {friendly name} {ff version}`.

## Special Cases and Exceptions

Use these to avoid repeated misunderstandings:

- **Force-field release DOI**:
  - Required by standards.
  - If PR notes DOI will be added before merge, keep as a single “must verify before merge” item.
  - Resolve once DOI is present.

- **Initial submission changelog**:
  - Optional for initial versions (`v4.0`, `v0.0`).
  - Do not flag presence of an optional changelog section as an issue.

- **Maintainer-accepted exceptions**:
  - If maintainers explicitly justify a non-default naming/wording choice as standards-compliant for that PR context, do not keep reopening it.

## Comment Style and Format

Use concise review output:
1. One-line experimental warning.
2. One short summary of compliance status.
3. Bullet list of findings with severity tags: `[blocking]`, `[non-blocking]`, `[manual-check]`.
4. Each finding must include:
   - path(s),
   - requirement being checked,
   - exact fix guidance.

If fully compliant, reply:
**“Looks good — all mandatory fields present and versioning correct. 🔥”**

## Comment Resolution Behavior

- Mark resolved when the file changed to address the issue.
- Mark resolved when a maintainer explanation demonstrates standards compliance.
- Do not keep unresolved comments that are based on misread requirements.
- Prefer one consolidated comment per issue category over many near-duplicates.