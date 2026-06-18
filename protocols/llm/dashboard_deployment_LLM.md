# LLM Protocol: Dashboard Serving & Static Asset Deployment

This protocol specifies the port mapping, background processes, log redirects, and static export procedures for the TOGA dashboards.

---

## 1. Running Dashboard Servers in the Background

Execute the backend python processes in the background (using nohup or nohup/nohup-like redirection):

### A. FastAPI TOGA MEME Dashboard (Port 8082)
* **Command**:
  ```bash
  python3 -u scripts/serve_meme_dashboard_fastapi.py 8082 > logs_serve_meme.log 2>&1 &
  ```
* **Lifespan/Startup Events**:
  1. Checks database index: `CREATE INDEX IF NOT EXISTS idx_site_gene_idx ON site_results(gene_name, site_index);`
  2. Runs background warm-up: queries `SELECT DISTINCT gene_name FROM site_results` and cache-warms key statistics.
  3. Listens on `0.0.0.0:8082`.

### B. TOGA aBSREL Dashboard (Port 8083)
* **Command**:
  ```bash
  python3 -u scripts/serve_absrel_dashboard.py 8083 > logs_serve_absrel.log 2>&1 &
  ```
* **Lifespan/Startup Events**:
  1. Precomputes global gene summaries.
  2. Listens on `localhost:8083`.

### C. TOGA Alignments Quality Dashboard (Port 8000)
* **Command**:
  ```bash
  python3 -u scripts/serve_alignments_dashboard.py 8000 > logs_serve_alignments.log 2>&1 &
  ```
* **Lifespan/Startup Events**:
  1. Runs robust outlier detection using log10-IQR & log10-Z-score on alignments database.
  2. Listens on `localhost:8000`.

---

## 2. Compiling and Exporting Static Assets
For static, database-free serving:

### A. MEME Static Site Generation
Run the compiler script:
```bash
python3 scripts/export_static_dashboard.py
```
* **Output Folder**: `toga-meme/`
* **Structure**:
  * `toga-meme/index.html`: Main dashboard layout with client-side table.
  * `toga-meme/data/summary.json`: Cached summary counts, metrics, and FDR thresholds.
  * `toga-meme/data/genes.json`: Complete list of gene metadata.
  * `toga-meme/data/alignments/{gene_name}.json`: Sliced sequence segments, codon grids, and site values.

### B. aBSREL Static Site Generation
Run the compiler script:
```bash
python3 scripts/export_absrel_static_dashboard.py
```
* **Output Folder**: `toga-absrel/`

### C. Synchronization to Server Share
Sync files to silverback web server host:
```bash
rsync -avz toga-meme/ silverback:/archive/sb-data/shares/web/web/toga-meme/
rsync -avz toga-absrel/ silverback:/archive/sb-data/shares/web/web/toga-absrel/
```
Ensure permission mask allows group/public read access on the target directory.
