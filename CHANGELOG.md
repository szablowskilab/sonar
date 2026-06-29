## 0.1.0.alpha.2

### New
- Add mapper module for mm2 based transcriptome mapping.
- New commands:
  - `map`: minimap2 based transcriptome mapping.
  - `append`: Appends forward and reverse adapter sequences to sensor sequences.
- Add minimap2 dependency.

### Breaking Changes

- Rename `generate_ses_lib` to `generate_lib`.
- Refactor binary and library internally to seperate crates.
- Rename sensors module to design and make it private.
  - Expose to public API from crate root.

## 0.1.0.alpha.1

- Initial library and CLI for generating RNA sensor libraries.
