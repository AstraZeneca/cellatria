## CellExpress Benchmark

This example demonstrates how to execute the CellExpress pipeline using the publicly available **GSE204716** dataset.

--- 

### Step 1: Download the dataset
Visit the GEO accession page for **GSE204716** and manually download the count matrices:
- [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204716](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204716)

Place the downloaded files into a local directory, e.g.:

```bash
mkdir -p ~/cellexpress_data/GSE204716
# Move your downloaded files here
```

--- 

### Step 2: Prepare the config file

Download the example configuration provided in this repository: [config_cellexpress.json](https://github.com/AstraZeneca/cellatria/blob/main/docs/config_cellexpress.json). 

Move it into your data directory (`~/cellexpress_data`)

--- 

### Step 3: Download required annotation models

Before running CellExpress, make sure to download the required pre-trained models for annotation:

**SCimilarity model (v1.1)**: Place it in: `~/scimilarity/model_v1.1`

> Read [SCimilarity Configuration](https://github.com/AstraZeneca/cellatria/blob/main/cellexpress/README.md#:~:text=Downloading%20Pretrained%20SCimilarity%20Models)

**CellTypist model (v1.6.3)**: Place it in: `~/celltypist/model_v1.6.3`

> Read [CellTypist Configuration](https://github.com/AstraZeneca/cellatria/blob/main/cellexpress/README.md#:~:text=Downloading%20Pretrained%20CellTypist%20Models)

**Note**: These paths must match the `input`, `sci_model_path`, and `cty_model_path` fields in your config file.
If needed, modify the config to reflect your local download paths.

--- 

### Step 4: Run CellExpress using Docker

Pull the container

```bash
docker pull ghcr.io/astrazeneca/cellatria:v1.0.0
```

Then, execute the CellExpress pipeline directly:

```bash
docker run -it --rm \
  -v ~/cellexpress_data/GSE204716:/data \
  ghcr.io/astrazeneca/cellatria:v1.0.0 cellexpress \
  --config ~/cellexpress_data/config_cellexpress_v1-0-0_ziFzM4S.json
```