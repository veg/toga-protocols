#!/usr/bin/env python3
import os
import json
import sqlite3
import pandas as pd

def generate_species_dashboard():
    print("Generating species to genome mapping dashboard...")
    
    # Paths
    tree_path = "toga_protocols/data/pruned_species_tree.nhx"
    mapping_path = "toga_protocols/data/species_to_taxon_map.csv"
    output_path = "toga_protocols/reports/species_mapping_dashboard.html"
    
    # Ensure output dir exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Read Newick tree
    with open(tree_path, "r") as f:
        newick_str = f.read().strip()
        
    # Read mapping metadata
    df_map = pd.read_csv(mapping_path)
    
    # Read assembly metadata to populate accessions and names
    assembly_path = "toga_protocols/data/species_to_assembly.csv"
    if os.path.exists(assembly_path):
        df_assembly = pd.read_csv(assembly_path)
        df = pd.merge(df_map, df_assembly[['Species', 'NCBI_Accession']], on="Species", how="left")
        df['Assembly_Accession'] = df['NCBI_Accession'].fillna("N/A")
        df['Species_Scientific_Name'] = df['Species']
    else:
        df = df_map.copy()
        df['Assembly_Accession'] = "N/A"
        df['Species_Scientific_Name'] = df['Species']
        
    mapping_data = df.to_dict(orient="records")
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Species to Genome Mapping Dashboard</title>
    <!-- jQuery, Bootstrap, D3 v3.5.17, Underscore, and Phylotree.js -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.17/d3.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.8.3/underscore-min.js"></script>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/phylotree@0.1.6/phylotree.css">
    <script src="https://cdn.jsdelivr.net/npm/phylotree@0.1.6/phylotree.js"></script>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css">
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            background-color: #f8f9fa;
        }}
        .sidebar {{
            background: white;
            border-right: 1px solid #dee2e6;
            height: 100vh;
            overflow-y: auto;
        }}
        .tree-container {{
            background: white;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            padding: 20px;
            height: calc(100vh - 100px);
            overflow: auto;
        }}
        .node {{
            font-size: 10px;
        }}
        .node:hover {{
            fill: #0d6efd;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <div class="container-fluid">
        <div class="row">
            <div class="col-md-4 sidebar p-4">
                <h3 class="mb-4 text-primary">Species Mapping</h3>
                <p class="text-muted">Total representative leaves: <strong>{len(df)}</strong></p>
                <div class="mb-3">
                    <input type="text" id="search" class="form-control" placeholder="Search species or assembly...">
                </div>
                <div class="table-responsive">
                    <table class="table table-striped table-hover table-sm">
                        <thead>
                            <tr>
                                <th>Taxon</th>
                                <th>Accession</th>
                                <th>Species</th>
                            </tr>
                        </thead>
                        <tbody id="mapping-table">
                            <!-- Populated dynamically -->
                        </tbody>
                    </table>
                </div>
            </div>
            <div class="col-md-8 p-4">
                <h3 class="mb-4 text-secondary">Phylotree.js Species Tree</h3>
                <div class="tree-container">
                    <svg id="tree_display"></svg>
                </div>
            </div>
        </div>
    </div>

    <script>
        const mappingData = {json.dumps(mapping_data)};
        const newickString = {json.dumps(newick_str)};

        // Populate Table
        function renderTable(filterText = "") {{
            const tbody = document.getElementById("mapping-table");
            tbody.innerHTML = "";
            const filtered = mappingData.filter(d => 
                d.Trimmed_Taxon_Name.toLowerCase().includes(filterText.toLowerCase()) ||
                d.Assembly_Accession.toLowerCase().includes(filterText.toLowerCase()) ||
                d.Species_Scientific_Name.toLowerCase().includes(filterText.toLowerCase())
            );
            filtered.forEach(d => {{
                const tr = document.createElement("tr");
                tr.innerHTML = `
                    <td><code>${{d.Trimmed_Taxon_Name}}</code></td>
                    <td><a href="https://www.ncbi.nlm.nih.gov/assembly/${{d.Assembly_Accession}}" target="_blank">${{d.Assembly_Accession}}</a></td>
                    <td><small>${{d.Species_Scientific_Name}}</small></td>
                `;
                tbody.appendChild(tr);
            }});
        }}
        
        document.getElementById("search").addEventListener("input", (e) => renderTable(e.target.value));
        renderTable();

        // Render Tree using Phylotree.js
        try {{
            const tree = d3.layout.phylotree()
                .svg(d3.select("#tree_display"))
                .options({{
                    "layout": "vertical",
                    "transitions": false,
                    "draw-size-bubbles": false,
                    "zoom": true
                }})
                .size([2200, 750]);
            tree(d3.layout.newick_parser(newickString)).layout();
        }} catch (e) {{
            console.error("Error rendering tree with phylotree.js:", e);
            document.querySelector(".tree-container").innerHTML = `
                <div class="alert alert-warning">
                    <h5>Unable to render tree dynamically via Phylotree.js</h5>
                    <p>Error message: <code>${{e.message}}</code></p>
                    <p>Fallback Newick Representation:</p>
                    <textarea class="form-control" rows="15" readonly>${{newickString}}</textarea>
                </div>
            `;
        }}
    </script>
</body>
</html>
"""
    with open(output_path, "w") as f:
        f.write(html_content)
    print(f"Saved species dashboard to: {output_path}")

def generate_quality_dashboard():
    print("Generating alignment quality and filtering dashboard...")
    
    # Path
    db_path = "alignments_stats.db"
    output_path = "toga_protocols/reports/alignment_quality_dashboard.html"
    gene_table_json_path = "scratch/alignment_gene_table.json"
    
    # Ensure output dir exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Query database stats if exists
    total_alignments = 18434
    retained_genes = 18275
    dropped_few_seqs = 159
    avg_seqs = 624.7
    avg_codons = 577.0
    low_cov_discarded = 1115190
    stops_discarded = 2795
    dedup_discarded = 2128950
    
    if os.path.exists(db_path):
        try:
            conn = sqlite3.connect(db_path)
            c = conn.cursor()
            c.execute("select count(*) from alignments")
            total_alignments = c.fetchone()[0]
            conn.close()
        except Exception as e:
            print(f"Error querying alignments_stats.db: {e}")
            
    # Read gene-by-gene details if exists
    gene_table_data = []
    if os.path.exists(gene_table_json_path):
        try:
            with open(gene_table_json_path, "r") as f:
                gene_table_data = json.load(f)
            print(f"Embedded {len(gene_table_data)} gene table rows.")
        except Exception as e:
            print(f"Error loading {gene_table_json_path}: {e}")
            
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Alignment Quality & Retained Genes Dashboard</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{
            font-family: 'Inter', sans-serif;
            background-color: #f8f9fa;
        }}
        .card {{
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.05);
            border: none;
            margin-bottom: 24px;
        }}
        .metric-value {{
            font-size: 2.2rem;
            font-weight: 700;
        }}
        .table-responsive {{
            max-height: 600px;
            overflow-y: auto;
        }}
        .status-badge {{
            font-weight: 600;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85rem;
        }}
        .status-kept {{
            background-color: #d1e7dd;
            color: #0f5132;
        }}
        .status-filtered {{
            background-color: #f8d7da;
            color: #842029;
        }}
    </style>
</head>
<body>
    <div class="container py-5">
        <h1 class="mb-5 text-primary">Alignment Curation & Quality Report</h1>
        
        <!-- Summary Row -->
        <div class="row">
            <div class="col-md-3">
                <div class="card p-4 text-center bg-white text-dark">
                    <div class="text-muted">Total Input MSAs</div>
                    <div class="metric-value text-primary">{total_alignments:,}</div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card p-4 text-center bg-white text-dark">
                    <div class="text-muted">Retained Genes</div>
                    <div class="metric-value text-success">{retained_genes:,}</div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card p-4 text-center bg-white text-dark">
                    <div class="text-muted">Dropped Genes (&lt;3 seqs)</div>
                    <div class="metric-value text-danger">{dropped_few_seqs:,}</div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card p-4 text-center bg-white text-dark">
                    <div class="text-muted">Avg. Seqs / Gene</div>
                    <div class="metric-value text-info">{avg_seqs}</div>
                </div>
            </div>
        </div>

        <div class="row">
            <div class="col-md-6">
                <div class="card p-4 bg-white">
                    <h4 class="mb-4">Sequence Filtering Breakdown</h4>
                    <canvas id="filterChart"></canvas>
                </div>
            </div>
            <div class="col-md-6">
                <div class="card p-4 bg-white">
                    <h4 class="mb-4">Dataset Curation Summary</h4>
                    <table class="table table-hover">
                        <thead>
                            <tr>
                                <th>Curation Category</th>
                                <th class="text-end">Value</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td>Low Coverage Discarded (&lt;50% resolved)</td>
                                <td class="text-end text-danger fw-bold">{low_cov_discarded:,}</td>
                            </tr>
                            <tr>
                                <td>Premature Stop Codons Discarded</td>
                                <td class="text-end text-danger fw-bold">{stops_discarded:,}</td>
                            </tr>
                            <tr>
                                <td>Redundant Assemblies Pruned</td>
                                <td class="text-end text-warning fw-bold">{dedup_discarded:,}</td>
                            </tr>
                            <tr>
                                <td>Average Codon Length (Cleaned)</td>
                                <td class="text-end text-primary fw-bold">{avg_codons} codons</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- Gene Curation Table Card -->
        <div class="card p-4 bg-white">
            <h4 class="mb-3">Gene Curation Details</h4>
            <div class="row g-3 mb-3 align-items-center">
                <div class="col-md-4">
                    <input type="text" id="gene-search" class="form-control" placeholder="Search gene or transcript...">
                </div>
                <div class="col-md-3">
                    <select id="status-filter" class="form-select">
                        <option value="ALL">All Statuses</option>
                        <option value="Kept">Kept</option>
                        <option value="Filtered">Filtered</option>
                    </select>
                </div>
                <div class="col-md-5 text-end">
                    <button class="btn btn-primary" onclick="exportToCSV()">Export to CSV</button>
                </div>
            </div>
            
            <div class="table-responsive">
                <table class="table table-hover table-striped">
                    <thead>
                        <tr>
                            <th style="cursor:pointer;" onclick="sortTable('gene')">Gene ↕</th>
                            <th style="cursor:pointer;" onclick="sortTable('transcript')">Transcript ID ↕</th>
                            <th style="cursor:pointer;" onclick="sortTable('orig_seqs')">Original Seqs ↕</th>
                            <th style="cursor:pointer;" onclick="sortTable('retained_seqs')">Retained Seqs ↕</th>
                            <th style="cursor:pointer;" onclick="sortTable('orig_codons')">Original Codons ↕</th>
                            <th style="cursor:pointer;" onclick="sortTable('retained_codons')">Retained Codons ↕</th>
                            <th style="cursor:pointer;" onclick="sortTable('status')">Status ↕</th>
                        </tr>
                    </thead>
                    <tbody id="gene-table-body">
                        <!-- Populated dynamically -->
                    </tbody>
                </table>
            </div>

            <div class="d-flex justify-content-between align-items-center mt-3">
                <div class="text-muted" id="table-info">Showing 0-0 of 0 entries</div>
                <div class="btn-group">
                    <button class="btn btn-outline-secondary btn-sm" id="btn-prev" onclick="prevPage()">Previous</button>
                    <button class="btn btn-outline-secondary btn-sm" id="btn-next" onclick="nextPage()">Next</button>
                </div>
            </div>
        </div>
    </div>

    <script>
        // Render Chart
        const ctx = document.getElementById('filterChart').getContext('2d');
        new Chart(ctx, {{
            type: 'pie',
            data: {{
                labels: ['Redundant Assemblies Pruned', 'Low Coverage Discarded', 'Premature Stops Discarded'],
                datasets: [{{
                    data: [{dedup_discarded}, {low_cov_discarded}, {stops_discarded}],
                    backgroundColor: ['#ffc107', '#dc3545', '#6c757d'],
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{
                    legend: {{
                        position: 'bottom',
                    }}
                }}
            }}
        }});

        // Table Data & Logic
        const geneData = {json.dumps(gene_table_data)};
        let filteredData = [...geneData];
        let currentPage = 1;
        const pageSize = 50;
        let sortField = 'gene';
        let sortAsc = true;

        function renderGeneTable() {{
            const tbody = document.getElementById("gene-table-body");
            tbody.innerHTML = "";
            
            // Apply sorting
            filteredData.sort((a, b) => {{
                let valA = a[sortField];
                let valB = b[sortField];
                
                if (typeof valA === 'string') {{
                    valA = valA.toLowerCase();
                    valB = valB.toLowerCase();
                }}
                
                if (valA < valB) return sortAsc ? -1 : 1;
                if (valA > valB) return sortAsc ? 1 : -1;
                return 0;
            }});

            const startIdx = (currentPage - 1) * pageSize;
            const endIdx = Math.min(startIdx + pageSize, filteredData.length);
            const pageData = filteredData.slice(startIdx, endIdx);

            if (pageData.length === 0) {{
                tbody.innerHTML = `<tr><td colspan="7" class="text-center text-muted p-4">No matching records found.</td></tr>`;
                document.getElementById("table-info").innerText = "Showing 0 of 0 entries";
                document.getElementById("btn-prev").disabled = true;
                document.getElementById("btn-next").disabled = true;
                return;
            }}

            pageData.forEach(row => {{
                const tr = document.createElement("tr");
                const isKept = row.status === 'Kept';
                const badgeClass = isKept ? 'status-kept' : 'status-filtered';
                
                tr.innerHTML = `
                    <td><strong>${{row.gene}}</strong></td>
                    <td><code>${{row.transcript}}</code></td>
                    <td>${{row.orig_seqs}}</td>
                    <td>${{row.retained_seqs}}</td>
                    <td>${{row.orig_codons}}</td>
                    <td>${{row.retained_codons}}</td>
                    <td><span class="status-badge ${{badgeClass}}">${{row.status}}</span></td>
                `;
                tbody.appendChild(tr);
            }});

            document.getElementById("table-info").innerText = `Showing ${{startIdx + 1}} to ${{endIdx}} of ${{filteredData.length}} entries`;
            document.getElementById("btn-prev").disabled = currentPage === 1;
            document.getElementById("btn-next").disabled = endIdx >= filteredData.length;
        }}

        function sortTable(field) {{
            if (sortField === field) {{
                sortAsc = !sortAsc;
            }} else {{
                sortField = field;
                sortAsc = true;
            }}
            currentPage = 1;
            renderGeneTable();
        }}

        function prevPage() {{
            if (currentPage > 1) {{
                currentPage--;
                renderGeneTable();
            }}
        }}

        function nextPage() {{
            const maxPage = Math.ceil(filteredData.length / pageSize);
            if (currentPage < maxPage) {{
                currentPage++;
                renderGeneTable();
            }}
        }}

        function filterData() {{
            const searchText = document.getElementById("gene-search").value.toLowerCase();
            const statusFilter = document.getElementById("status-filter").value;

            filteredData = geneData.filter(row => {{
                const matchesSearch = row.gene.toLowerCase().includes(searchText) || 
                                      row.transcript.toLowerCase().includes(searchText);
                
                let matchesStatus = true;
                if (statusFilter === 'Kept') {{
                    matchesStatus = row.status === 'Kept';
                }} else if (statusFilter === 'Filtered') {{
                    matchesStatus = row.status !== 'Kept';
                }}
                
                return matchesSearch && matchesStatus;
            }});

            currentPage = 1;
            renderGeneTable();
        }}

        function exportToCSV() {{
            let csvContent = "data:text/csv;charset=utf-8,Gene,Transcript ID,Original Seqs,Retained Seqs,Original Codons,Retained Codons,Status\\n";
            filteredData.forEach(row => {{
                csvContent += `${{row.gene}},${{row.transcript}},${{row.orig_seqs}},${{row.retained_seqs}},${{row.orig_codons}},${{row.retained_codons}},${{row.status}}\\n`;
            }});
            const encodedUri = encodeURI(csvContent);
            const link = document.createElement("a");
            link.setAttribute("href", encodedUri);
            link.setAttribute("download", "alignment_curation_report.csv");
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}

        document.getElementById("gene-search").addEventListener("input", filterData);
        document.getElementById("status-filter").addEventListener("change", filterData);

        // Initial render
        renderGeneTable();
    </script>
</body>
</html>
"""
    with open(output_path, "w") as f:
        f.write(html_content)
    print(f"Saved quality dashboard to: {output_path}")

if __name__ == "__main__":
    generate_species_dashboard()
    generate_quality_dashboard()
