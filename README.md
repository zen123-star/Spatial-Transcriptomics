# Spatial-Transcriptomics
# Spatial Transcriptomics Analysis

#*FILES:*
1.*Spatial Transcriptomics*
https://drive.google.com/file/d/1FwM2L8RhW5pUo1hbW9gOUSR9YNhP6GV3/view?usp=sharing

2.*Visium flouresecnce*
https://drive.google.com/file/d/10vppxJ2M_F5DueE3EFm-vvOaGG7LJr2j/view?usp=sharing

3.*Visisum H&E data:*
https://drive.google.com/file/d/10vppxJ2M_F5DueE3EFm-vvOaGG7LJr2j/view?usp=sharing

4.*Xenium data*
https://drive.google.com/file/d/10vppxJ2M_F5DueE3EFm-vvOaGG7LJr2j/view?usp=sharing

A comprehensive collection of four end-to-end tutorials covering spatial transcriptomics data analysis using **Scanpy**, **Squidpy**, and **SpatialData**. Each tutorial targets a distinct platform or data modality, from classical Visium H&E slides to next-generation Xenium single-cell resolved data.

---

## Table of Contents

1. [Tutorial 1 — Scanpy: Basic Visium & MERFISH Analysis](#tutorial-1--scanpy-basic-visium--merfish-analysis)
2. [Tutorial 2 — Squidpy: Visium Fluorescence Image Analysis](#tutorial-2--squidpy-visium-fluorescence-image-analysis)
3. [Tutorial 3 — Squidpy: Visium H&E Data with Spatial Statistics](#tutorial-3--squidpy-visium-he-data-with-spatial-statistics)
4. [Tutorial 4 — Squidpy + SpatialData: Xenium Single-Cell Analysis](#tutorial-4--squidpy--spatialdata-xenium-single-cell-analysis)

---

## Tutorial 1 — Scanpy: Basic Visium & MERFISH Analysis

### Overview

This tutorial provides a foundational introduction to spatial transcriptomics analysis using **Scanpy**. It demonstrates the full pipeline on a **10x Genomics Visium** dataset of a human lymph node and includes a secondary example using **MERFISH** data from Xia et al. (2019).

### Dataset

| Property | Details |
|---|---|
| **Platform** | 10x Genomics Visium |
| **Tissue** | Human Lymph Node |
| **Cells / Spots** | 4,035 spots × 36,601 genes (pre-filter) |
| **Secondary Dataset** | MERFISH — cultured U2-OS cells (645 cells × 12,903 genes) |

### Workflow

#### 1. Data Loading & QC
- Load Visium data via `sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")`
- Calculate QC metrics: total counts, genes per spot, mitochondrial read percentage using `sc.pp.calculate_qc_metrics()`
- Visualize count distributions with histograms to determine filtering thresholds

#### 2. Preprocessing & Filtering
- Filter spots by count range: `min_counts=5000`, `max_counts=35000`
- Exclude spots with >20% mitochondrial counts
- Filter low-coverage genes: `min_cells=10`
- Normalize counts per cell and apply log1p transformation
- Identify top 2,000 highly variable genes using the Seurat flavor

#### 3. Dimensionality Reduction & Clustering
- PCA (50 components) → KNN graph construction → UMAP embedding
- Leiden clustering (igraph flavor, undirected, 2 iterations) yields **10 clusters**

#### 4. Spatial Visualization
- Overlay cluster labels on H&E tissue image using `sc.pl.spatial()` with the `hires` image key
- Zoom into regions of interest using `crop_coord` to examine clusters 5 and 9
- Adjust spot transparency (`alpha`) to reveal underlying tissue morphology

  <img width="2038" height="609" alt="image" src="https://github.com/user-attachments/assets/df464b0f-eae8-4403-85db-18b9bcc2cff6" />

<img width="2294" height="1091" alt="image" src="https://github.com/user-attachments/assets/afe6eeab-9450-4a1b-bcfa-4ef6848d3471" />


#### 5. Marker Gene Analysis
- Rank marker genes per cluster using `sc.tl.rank_genes_groups()` (t-test method)
- Visualize top-10 markers for cluster 9 via a heatmap
- Spatially map key marker genes (e.g., `CR2`, `COL1A2`, `SYPL1`) onto the tissue

#### 6. MERFISH Extension
- Load coordinate and count matrices from the original publication
- Assign spatial coordinates to `adata.obsm["spatial"]`
- Perform standard preprocessing, PCA (15 components), and Leiden clustering (resolution=0.5) → **6 clusters** (corresponding to cell-cycle stages)
- Visualize clusters in UMAP and spatial coordinate space using `sc.pl.embedding()`

### Key Functions

```python
sc.datasets.visium_sge()         # Load Visium dataset
sc.pp.calculate_qc_metrics()     # QC metrics computation
sc.pp.normalize_total()          # Count normalization
sc.pp.highly_variable_genes()    # HVG selection
sc.tl.leiden()                   # Graph-based clustering
sc.pl.spatial()                  # Spatial overlay visualization
sc.tl.rank_genes_groups()        # Marker gene identification
sc.pl.embedding(basis="spatial") # MERFISH spatial plot
```

### Dependencies

```
scanpy >= 1.10.0
anndata >= 0.11.0
matplotlib, seaborn, pandas
```

---

## Tutorial 2 — Squidpy: Visium Fluorescence Image Analysis

### Overview

This tutorial demonstrates the use of **Squidpy's image analysis capabilities** on a **Visium fluorescence** dataset. The focus is on extracting multi-scale image features from the tissue image and deriving morphology-based cluster annotations that complement gene expression clusters. The dataset is a cropped coronal section of the **mouse brain**.

### Dataset

| Property | Details |
|---|---|
| **Platform** | 10x Genomics Visium (Fluorescence) |
| **Tissue** | Mouse Brain (coronal section, cropped) |
| **Image Channels** | DAPI (DNA), anti-NEUN (neurons), anti-GFAP (glial cells) |
| **Spots** | 704 Visium spots |

### Workflow

#### 1. Data Loading
- Load pre-processed AnnData and ImageContainer using `sq.datasets.visium_fluo_adata_crop()` and `sq.datasets.visium_fluo_image_crop()`
- Visualize pre-annotated gene-expression clusters with `sq.pl.spatial_scatter()`
- Inspect fluorescence channels with `img.show(channelwise=True)`
<img width="790" height="289" alt="image" src="https://github.com/user-attachments/assets/c6701bed-99e6-491d-9ff8-1f388905ccaf" />

#### 2. Image Segmentation
- Smooth the image layer using `sq.im.process(method="smooth")`
- Perform nucleus segmentation on the DAPI channel (channel 0) via watershed algorithm: `sq.im.segment(method="watershed")`
- Results stored in `img["segmented_watershed"]` as a label image
<img width="515" height="267" alt="image" src="https://github.com/user-attachments/assets/41c7310a-38b6-4554-9121-093ac21960af" />

#### 3. Segmentation Feature Extraction
- Extract segmentation features with `sq.im.calculate_image_features(features="segmentation")`
- Features include:
  - **Cell count per spot** — approximation of cellular density
  - **Channel-wise mean intensity** — per-channel signal within segmented objects
- Use `sq.pl.extract()` to temporarily project `obsm` features into `obs` for plotting

<img width="1088" height="829" alt="image" src="https://github.com/user-attachments/assets/2aa0d019-e197-4dc7-a5be-63fc72d9eb56" />

#### 4. Multi-Scale Image Feature Extraction
Three feature sets extracted at different scales:

| Feature Set | Features Included | Scale | Context |
|---|---|---|---|
| `features_orig` | summary, texture, histogram | 1.0 | Spot only (masked circle) |
| `features_context` | summary, histogram | 1.0 | Full spot crop |
| `features_lowres` | summary, histogram | 0.25 | Low-resolution context |

All feature sets are concatenated into `adata.obsm["features"]` with deduplicated column names.

#### 5. Feature-Based Clustering
- A helper `cluster_features()` function applies PCA → KNN → Leiden clustering on any feature subset
- Generate three cluster annotations:
  - `features_summary_cluster`
  - `features_histogram_cluster`
  - `features_texture_cluster`
- Compare image-feature clusters against gene-expression clusters spatially
<img width="3435" height="2205" alt="image" src="https://github.com/user-attachments/assets/7d982724-43d2-4462-b5ea-16ca248d6b8b" />

### Key Findings

- **Cell-rich pyramidal layer** of the Hippocampus is resolved at higher granularity by image features than by gene expression alone
- **anti-NEUN intensity** (channel 1) highlights Cortex_1 and Cortex_3 as neuron-enriched regions
- **anti-GFAP intensity** (channel 2) marks Fiber_tracts and lateral ventricles as glial-rich areas

### Key Functions

```python
sq.im.process()                      # Image smoothing / preprocessing
sq.im.segment()                      # Cell/nucleus segmentation
sq.im.calculate_image_features()     # Feature extraction pipeline
sq.pl.extract()                      # Project obsm features into obs
sq.pl.spatial_scatter()              # Spatial visualization
```

### Dependencies

```
squidpy >= 1.2.3
scanpy, anndata
pandas, matplotlib
```

---

## Tutorial 3 — Squidpy: Visium H&E Data with Spatial Statistics

### Overview

This tutorial applies **Squidpy's spatial statistics and graph analysis tools** to a Visium H&E dataset. Starting from pre-annotated gene-expression clusters, it demonstrates neighborhood enrichment, co-occurrence scoring, ligand-receptor interaction analysis, and spatially variable gene detection. The dataset is a coronal section of the **mouse brain**.

### Dataset

| Property | Details |
|---|---|
| **Platform** | 10x Genomics Visium (H&E stained) |
| **Tissue** | Mouse Brain (coronal section) |
| **Spots** | 2,688 Visium spots |
| **Annotations** | Pre-annotated clusters (Allen Brain Atlas, Linnarsson lab) |

### Workflow

#### 1. Data Loading & Image Feature Extraction
- Load data via `sq.datasets.visium_hne_adata()` and `sq.datasets.visium_hne_image()`
- Extract multi-scale summary features at two scales (`scale=1.0` and `scale=2.0`) covering increasing spatial context
- Concatenate into `adata.obsm["features"]` and compute feature-based Leiden clusters
- Compare image-feature clusters vs. gene-expression clusters spatially
<img width="651" height="221" alt="image" src="https://github.com/user-attachments/assets/50671b8c-5762-48c0-b0a5-86d4350d0a39" />

#### 2. Spatial Graph Construction
- Build a spatial connectivity matrix using `sq.gr.spatial_neighbors()`
- Required for all downstream graph-based statistics

#### 3. Neighborhood Enrichment Analysis
- Compute permutation-based enrichment scores with `sq.gr.nhood_enrichment(n_perms=1000)`
- Score reflects how often spots from two clusters are spatial neighbors vs. random expectation
- **Result:** High enrichment between `Pyramidal_layer`, `Pyramidal_layer_dentate_gyrus`, and the broader `Hippocampus` cluster

#### 4. Co-Occurrence Analysis
- Compute conditional probability ratios across increasing spatial radii using `sq.gr.co_occurrence()`
- Formula: `p(exp | cond) / p(exp)` — score > 1 indicates spatial enrichment
- **Result:** Pyramidal_layer co-occurs at short distances with the Hippocampus cluster

#### 5. Ligand-Receptor Interaction Analysis
- Re-implementation of **CellPhoneDB** extended with the **Omnipath** database
- Run with `sq.gr.ligrec(n_perms=100, cluster_key="cluster")`
- Visualize candidate interactions between `Hippocampus` (source) and `Pyramidal_layer` / `Pyramidal_layer_dentate_gyrus` (targets)
- Filter with `means_range=(3, np.inf)` and `alpha=1e-4` to select high-confidence interactions

#### 6. Spatially Variable Genes — Moran's I
- Evaluate spatial autocorrelation using `sq.gr.spatial_autocorr(mode="moran")`
- Applied to top 1,000 highly variable genes with 100 permutations
- Results saved to `adata.uns["moranI"]` sorted by Moran's I statistic

- <img width="2586" height="431" alt="image" src="https://github.com/user-attachments/assets/d6aebd27-c754-4d2b-a281-7451bc68254a" />


**Top spatially variable genes identified:**

| Gene | Moran's I |
|---|---|
| Olfm1 | 0.763 |
| Plp1 | 0.748 |
| Itpka | 0.727 |
| Snap25 | 0.721 |
| Nnat | 0.709 |

### Key Functions

```python
sq.gr.spatial_neighbors()       # Spatial connectivity graph
sq.gr.nhood_enrichment()        # Neighborhood enrichment score
sq.gr.co_occurrence()           # Spatial co-occurrence probability
sq.gr.ligrec()                  # Ligand-receptor interaction analysis
sq.gr.spatial_autocorr()        # Moran's I / Geary's C
sq.pl.nhood_enrichment()        # Enrichment heatmap
sq.pl.co_occurrence()           # Co-occurrence line plots
sq.pl.ligrec()                  # Dotplot of LR interactions
```

### Dependencies

```
squidpy >= 1.2.3
scanpy, anndata, numpy
pandas, matplotlib
```

---

## Tutorial 4 — Squidpy + SpatialData: Xenium Single-Cell Analysis

### Overview

This tutorial demonstrates analysis of **10x Genomics Xenium** data — a subcellular-resolution, in-situ sequencing platform. It uses the **SpatialData** framework for data ingestion and storage, **Scanpy** for transcriptomic analysis, and **Squidpy** for spatial statistics. The dataset is a **Human Lung Cancer** sample with 161,000 cells and 480 genes.

### Dataset

| Property | Details |
|---|---|
| **Platform** | 10x Genomics Xenium |
| **Tissue** | Human Lung Cancer |
| **Cells** | 161,000 cells × 480 genes |
| **Transcripts** | 40,257,199 (3D points) |
| **Image** | `morphology_focus` (multichannel, multiscale) |

### Workflow

#### 1. Data Ingestion with SpatialData
- Parse the raw Xenium output directory using the `xenium()` reader from `spatialdata-io`
- Convert and persist the full dataset to Zarr format for efficient I/O: `sdata.write(zarr_path)`
- The SpatialData object contains:
  - **Images:** `morphology_focus` (5 resolution levels)
  - **Labels:** `cell_labels`, `nucleus_labels`
  - **Points:** `transcripts` DataFrame
  - **Shapes:** `cell_boundaries`, `nucleus_boundaries`, `cell_circles`
  - **Tables:** `table` (AnnData with count matrix)

#### 2. Quality Control
- Compute QC metrics with `sc.pp.calculate_qc_metrics(percent_top=(10, 20, 50, 150))`
- Calculate negative control fractions:
  - Negative DNA probe count: ~0.005%
  - Negative decoding count: ~0.003%
- Visualize distributions of total transcripts, unique transcripts, cell area, and nucleus-to-cell area ratio

#### 3. Preprocessing & Embedding
- Filter cells: `min_counts=10`; filter genes: `min_cells=5`
- Preserve raw counts in `adata.layers["counts"]`
- Normalize → log1p → PCA → KNN neighbors → UMAP → Leiden clustering
- Visualize UMAP colored by total counts, genes-by-counts, and Leiden clusters
- Visualize Leiden clusters in spatial coordinates with `sq.pl.spatial_scatter(shape=None)`

#### 4. Centrality Scores
- Build spatial graph: `sq.gr.spatial_neighbors(coord_type="generic", delaunay=True)`
- Compute three centrality metrics per cluster with `sq.gr.centrality_scores()`:
  - **Closeness centrality** — how close a group is to all other nodes
  - **Degree centrality** — fraction of non-group nodes connected to the group
  - **Clustering coefficient** — degree to which group nodes cluster together
<img width="1561" height="431" alt="image" src="https://github.com/user-attachments/assets/02391258-f81b-48ac-9ddf-c04d2348ca01" />

#### 5. Co-Occurrence Probability
- Subsample 50% of cells for computational efficiency: `sc.pp.subsample(fraction=0.5)`
- Compute co-occurrence scores with `sq.gr.co_occurrence(cluster_key="leiden")`
- Visualize probability ratios across spatial radii for Leiden cluster "12"

#### 6. Neighborhood Enrichment
- Compute permutation-based enrichment (1,000 permutations): `sq.gr.nhood_enrichment()`
- Side-by-side visualization of enrichment heatmap and spatial scatter plot

<img width="1224" height="494" alt="image" src="https://github.com/user-attachments/assets/c6748ce4-1a94-4bac-a555-fb5ee2d7aa2b" />

#### 7. Moran's I — Spatially Variable Genes
- Rebuild spatial graph on subsampled data with Delaunay triangulation
- Run `sq.gr.spatial_autocorr(mode="moran", n_perms=100)` on all 480 genes

**Top spatially variable genes identified:**

| Gene | Moran's I | Biological Relevance |
|---|---|---|
| AREG | 0.696 | Epithelial growth factor |
| MET | 0.683 | Receptor tyrosine kinase |
| ANXA1 | 0.667 | Anti-inflammatory mediator |
| EPCAM | 0.633 | Epithelial cell marker |
| DMBT1 | 0.588 | Tumor suppressor |

#### 8. Visualization with SpatialData-Plot & napari
- Use `spatialdata_plot` to overlay gene expression on the morphology image:
  ```python
  sdata.pl.render_images("morphology_focus")
      .pl.render_shapes("cell_circles", color="AREG", table_name="table")
      .pl.show(coordinate_systems="global")
  ```
- Interactive exploration available via `napari-spatialdata`: `Interactive(sdata)`

### Key Functions

```python
xenium()                             # Parse Xenium raw output
sdata.write()                        # Persist to Zarr
sq.gr.spatial_neighbors()           # Delaunay spatial graph
sq.gr.centrality_scores()           # Closeness / degree / clustering
sq.gr.co_occurrence()               # Conditional co-occurrence score
sq.gr.nhood_enrichment()            # Permutation enrichment
sq.gr.spatial_autocorr()            # Moran's I autocorrelation
sq.pl.centrality_scores()           # Centrality visualization
sdata.pl.render_images()            # SpatialData-plot overlay
napari_spatialdata.Interactive()    # Interactive GUI
```

### Dependencies

```
squidpy >= 1.2.3
spatialdata, spatialdata-io, spatialdata-plot
scanpy, anndata
napari, napari-spatialdata
seaborn, matplotlib
```

---

## Environment Setup

### Conda (Recommended)

```bash
conda env create -f environment.yml
conda activate spatial-transcriptomics
```

### pip

```bash
pip install scanpy squidpy spatialdata spatialdata-io spatialdata-plot
pip install napari napari-spatialdata
pip install openpyxl  # Required for MERFISH Excel coordinate files
```

---

## Repository Structure

```
.
├── 01_scanpy_visium_merfish.ipynb       # Tutorial 1 — Scanpy Visium & MERFISH
├── 02_squidpy_visium_fluorescence.ipynb # Tutorial 2 — Squidpy Fluorescence
├── 03_squidpy_visium_hne.ipynb          # Tutorial 3 — Squidpy H&E & Spatial Stats
├── 04_squidpy_xenium.ipynb              # Tutorial 4 — Xenium + SpatialData
├── data/
│   ├── V1_Human_Lymph_Node/             # Auto-downloaded by scanpy
│   ├── pnas.1912459116.sd12.csv         # MERFISH counts
│   ├── pnas.1912459116.sd15.xlsx        # MERFISH coordinates
│   ├── Xenium/                          # Raw Xenium output directory
│   └── Xenium.zarr/                     # Converted Zarr store
├── environment.yml
└── README.md
```

---

## References

- **Scanpy:** Wolf et al. (2018). *Genome Biology*
- **Squidpy:** Palla et al. (2022). *Nature Methods*
- **SpatialData:** Marconato et al. (2024). *Nature Methods*
- **CellPhoneDB:** Efremova et al. (2020). *Nature Protocols*
- **Omnipath:** Türei et al. (2016). *Nature Methods*
- **MERFISH dataset:** Xia et al. (2019). *PNAS*
- **Allen Brain Atlas:** https://mouse.brain-map.org
- **10x Genomics Visium datasets:** https://support.10xgenomics.com/spatial-gene-expression/datasets

---

## Authors & Acknowledgements

Tutorial content authored by **Giovanni Palla** and the **scverse** community. Analysis pipelines reference methods from the Theis Lab and the Linnarsson Lab. Cluster annotations for mouse brain datasets were performed using the Allen Brain Atlas and the Mouse Brain Gene Expression Atlas.

---

*© 2026 scverse / Alex Wolf, Fidel Ramirez, Sergei Rybakov. Licensed under the Read the Docs stable build.*
