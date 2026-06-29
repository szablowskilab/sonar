## 0.1.0.alpha.2

### New
- Add specificity module for mm2 based transcriptome mapping.
- Add spec CLI subcommand.
- Add minimap2 dependency.

### Breaking Changes

- Rename `generate_ses_lib` to `generate_lib`.
- Refactor binary and library internally to seperate crates.
- Rename sensors module to design and make it private.
  - Expose to public API from crate root.

## 0.1.0.alpha.1

- Initial library and CLI for generating RNA sensor libraries.
