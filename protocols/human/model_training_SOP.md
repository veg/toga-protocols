# Standard Operating Procedure (SOP): AxoMeme Model Training

This standard operating procedure outlines the high-level workflow for training the Phylogenetic Axial Transformer (AxoMeme) regression model. This model is trained to predict evolutionary selection pressure directly from Multiple Sequence Alignments (MSAs).

---

## 1. Overview of AxoMeme Training
The AxoMeme model is a phylogenetic neural network that predicts continuous selection pressure values (Likelihood Ratio Test statistics) using amino acid alignment slices and phylogenetic branch distances. Training requires three core inputs:
* **MEME Results Database**: The primary dataset of codon site selection scores.
* **Curated MSA Directory**: gzipped NEXUS files containing trimmed codon alignments and phylogenetic trees.
* **Precomputed Distance Cache**: A serialized index containing patristic distance matrices and coordinate mappings.

---

## 2. Training Sequence and Execution Steps

### Step 1: Precomputing Feature Coordinates
Calculating phylogenetic branch lengths and species distances on the fly is computationally bottlenecked. Features must be precompiled:
1. Run the cache precomputation script (`precompute_msa_cache.py`) to process tree branch lengths.
2. The script extracts pairwise distances, computes Multi-Dimensional Scaling (MDS) coordinates, and flags amino-acid variable positions.

### Step 2: Training Run Execution
Execute training using the main regression script (`train_regression.py`):
1. **Gene-based Partitioning**: Split data into 80% training and 20% validation by gene. Grouping by gene prevents model memorization of homologous sites (homology leakage).
2. **Resampling**: Downsample neutral sites (with zero selection signal) to balance the classification target.
3. **Loss Function**: Train using a composite loss combining **Weighted Huber Loss** (to focus on selection hotspots) and **Batch Ranking Loss** (to optimize relative score rankings).

### Step 3: Checkpoint Selection
Compare validation metrics at the end of each training epoch:
* **Spearman Rank Correlation ($\rho$)**: The primary optimization metric measuring the alignment of predicted vs. MLE orderings.
* **Pearson Correlation ($R$)**: Linear correlation.
* **Mean Squared Error (MSE)**: Absolute regression error.

Save the model configuration with the highest validation Spearman correlation to `MEME_transformer_joint.pt`.

---

## 3. Deployment and Hardware Warnings

* **Apple Silicon Precision Bug**: Do not run AxoMeme inference on Apple Silicon GPU accelerators (`device="mps"`). Metal Performance Shaders (MPS) have known numerical precision bugs with deep transformer graphs that can distort selection predictions. Always run inference on CPU on macOS devices.
