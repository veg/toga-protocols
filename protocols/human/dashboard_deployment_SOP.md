# Standard Operating Procedure (SOP): Dashboard Serving & Deployment

This standard operating procedure describes the workflows for launching, caching, and exporting the three interactive dashboards: MEME selection dashboard, aBSREL selection dashboard, and Alignments quality dashboard.

---

## 1. Local Interactive Dashboard Servers
The project uses three local services to inspect genomic analyses:

1. **MEME Selection Dashboard (Port 8082)**:
   * Served via a FastAPI backend script (`serve_meme_dashboard_fastapi.py`).
   * Automatically optimizes indexes, precomputes summary stats, and warms up the species cache upon startup.
2. **aBSREL Selection Dashboard (Port 8083)**:
   * Served via a Python script (`serve_absrel_dashboard.py`).
   * Pre-calculates gene-level branch selection statistics and caches them in memory.
3. **Alignments Quality Dashboard (Port 8000)**:
   * Served via a Python script (`serve_alignments_dashboard.py`).
   * Evaluates sequence curation metrics and caches the outlier detection results.

---

## 2. Deploying Dashboard Servers in the Background
To start the dashboard servers and ensure they persist in the background:
1. Run each server redirecting output logs to their respective log files.
2. Verify they are successfully listening on their targeted ports (`8082`, `8083`, `8000`).

---

## 3. Exporting Static Dashboards for Public Hosting
For public sharing or presentation, compile the dynamic database results into static JSON assets and self-contained static HTML pages:

1. **Export the MEME Static Dashboard**:
   * Run the export script (`export_static_dashboard.py`).
   * This generates static index pages and alignment segment assets in the `toga-meme/` public directory.
2. **Export the aBSREL Static Dashboard**:
   * Run the export script (`export_absrel_static_dashboard.py`).
   * This generates static index pages and branch selection assets in the `toga-absrel/` public directory.
3. **Sync to Web Server**:
   * Synchronize the public-facing static folders (`toga-meme/` and `toga-absrel/`) to the web hosting directory on `silverback` using `rsync`.
   * The static dashboards are now publicly accessible.
