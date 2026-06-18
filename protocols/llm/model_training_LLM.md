# LLM Protocol: PhyloAxialTransformer (AxoMeme) Model Training

This protocol defines the algorithmic implementation details, training mechanics, and scripts for training the **PhyloAxialTransformer** (AxoMeme) regression model.

---

## 1. Feature Extraction & Coordinate Precomputation
Active coordinate calculation during GPU training is CPU-bound. Alignments and trees are compiled into an MDS-coordinate cache.

### A. Algorithmic Steps in `precompute_msa_cache.py`
1. **Tree Parsing**: Read the embedded Newick tree from each gzipped NEXUS file in `msa/`.
2. **Patristic Distance Computation**: For all species pairs $(i, j)$ in the tree, calculate the sum of branch lengths connecting them, yielding the distance matrix $D$.
3. **Multi-Dimensional Scaling (MDS)**: Apply classical MDS to reduce $D$ into a 4-dimensional spatial coordinate matrix $X \in \mathbb{R}^{S \times 4}$, where $S$ is the number of species (NTAX):
   $$B = -\frac{1}{2} H D^2 H \quad \text{where} \quad H = I - \frac{1}{n}\mathbf{1}\mathbf{1}^T$$
   Compute the eigendecomposition of $B$ and extract coordinates using the top 4 eigenvectors scaled by their eigenvalues.
4. **Variability Filter**: Identify codon sites showing nucleotide variability, excluding synonymous transitions in Serine codon islands (TCN vs. AGY).

### B. Execution
```bash
python3 scripts/precompute_msa_cache.py \
  --db_path /Users/sergei/Projects/TOGA_MEME/meme_results.db \
  --msa_dir /Users/sergei/Projects/TOGA_MEME/msa \
  --out_cache /Users/sergei/Projects/TOGA_MEME/scratch/msa_cache.pkl.gz
```

---

## 2. Model Training Protocol
Train the transformer model to predict selection pressure: $y = \ln(\max(0, \text{LRT}) + 1)$ using [train_regression.py](file:///Users/sergei/Dropbox/TOGA2026/scripts/train_regression.py).

### A. Algorithmic Configurations
* **Grouped Data Splitting**: Partition data 80% train, 20% validation. Ensure all codon sites belonging to the same gene are grouped entirely within either train or validation to prevent homology leakage.
* **Variability Filtering**: Discard completely invariant sites since their target LRT value is $0.0$.
* **Balanced Resampling**: In each epoch, downsample the zero-signal sites ($y = 0$) to match the total count of positive-signal sites ($y > 0$).
* **Composite Optimization Loss**:
  1. **Weighted Huber Loss**: Regression loss scaled by target magnitude to penalize errors on selection hotspots:
     $$\mathcal{L}_{huber} = \text{Huber}(y_{pred}, y_{true}) \cdot (y_{true} + 1.0)$$
  2. **Batch Ranking Loss (Margin)**: Optimizes relative order alignment:
     $$\mathcal{L}_{rank} = \sum_{i,j} \max(0, - \text{sign}(y_i - y_j)(p_i - p_j) + \text{margin})$$
     Optimize using SGD or AdamW with $\mathcal{L}_{total} = \mathcal{L}_{huber} + \alpha \mathcal{L}_{rank}$.

### B. Execution
```bash
python3 scripts/train_regression.py \
  --db_path /Users/sergei/Projects/TOGA_MEME/meme_results.db \
  --msa_dir /Users/sergei/Projects/TOGA_MEME/msa \
  --cache_path /Users/sergei/Projects/TOGA_MEME/scratch/msa_cache.pkl.gz \
  --epochs 10
```

---

## 3. Post-Training Downstream Inference

Use [predict_regression_nexus.py](file:///Users/sergei/Dropbox/TOGA2026/scripts/predict_regression_nexus.py) to execute inference on a target alignment:

```bash
python3 scripts/predict_regression_nexus.py \
  --alignment /path/to/target.nex \
  --model /Users/sergei/Projects/TOGA_MEME/MEME_transformer_joint.pt \
  --output /path/to/predictions.csv \
  --device cpu
```

> [!CAUTION]
> **Apple Silicon MPS Backend Limitation**: Do not use the MPS accelerator backend (`--device mps`) for inference. MPS kernels have known precision bugs with deep transformer graphs, leading to floating point errors in selection probability estimation. Always enforce CPU inference on macOS.
