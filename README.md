# Multimodal single-cell inference reveals cis-regulatory network remodeling

Analysis and figure code for the Scientific Reports paper.

## Structure

| Folder | Contents |
| --- | --- |
| `figure1/` | Result PDF only (analysis-flow schematic; plotting source not available) |
| `figure2/` – `figure6/` | Standalone script `figureN.py` + result PDF |
| `S1/`, `S2/`, `S4/` | Standalone script `SN.py` + result PDF |
| `S3/` | Result PDF only (partial plotting source missing) |

Each figure folder is self-contained: the script does not import or read any other
script in this repository. Originally shared helper modules have been inlined with
namespace prefixes (e.g. `network_`, `validation_`) so every script runs on its own.

## Requirements

See `requirements.txt`. Core: numpy, pandas, matplotlib, scipy, scikit-learn,
networkx, pillow; the single-cell scripts additionally require anndata and scanpy.

## Usage

Each script reads input data from `DATA_DIR` and writes results to `OUTPUT_ROOT`,
both set at the top of the file. They can also be overridden with the
`CRC_DATA_DIR` and `CRC_OUTPUT_DIR` environment variables.

```bash
python figure2.py
```

## Notes

- `figure1/` — the source for the analysis-flow schematic is not available; only
  the result PDF is provided.
- `S3/` — parts of the original plotting source are missing; only the result PDF
  is provided.
- Panels whose input data cannot be found are reported as missing in the script's
  `运行报告.json` rather than being replaced with placeholder figures.

## License

MIT. See `LICENSE`.
