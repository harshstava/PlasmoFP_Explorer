# PlasmoFP Explorer

A Streamlit web application for exploring PlasmoFP gene predictions and annotations across 19 *Plasmodium* species. For use with the associated manuscript - PlasmoFP: leveraging deep learning to predict protein function of uncharacterized proteins across the malaria parasite genus 

### Prerequisites
- Python 3.8+

### Installation

1. **Clone or navigate to the project directory:**
   ```bash
   cd "PlasmoFP_Explorer"
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

### Running the Application

```bash
streamlit run plasmoFP_explorer_simple.py
```

The application will open in your browser at `http://localhost:8501`

## Data Structure

The application reads PlasmoFP dictionaries from `with_PFP_predictions_complete_2/` containing:

- **19 *Plasmodium* species** gene dictionaries
- **PFP predictions** at 5 FDR thresholds (0.01, 0.05, 0.10, 0.20, 0.30)
- **Original GO annotations** (experimental + IEA)
- **TM-Vec protein embeddings** for external usage (not utilized here)
- **Gene ontology aspects**: Molecular Function (MF), Biological Process (BP), Cellular Component (CC)

## How to use

### Search Features

#### 1. Gene ID Search
- Enter a gene ID (e.g., `PF3D7_0001`, `PVX_0001`)
- Use partial matches (e.g., `PF3D7` to find all P. falciparum genes)
- Support for comma-separated multiple gene searches
- View detailed results in expandable sections

#### 2. Gene Product Descritipn
- Search by protein product descriptions
- Find genes with specific functional keywords
- Fuzzy matching for flexible searches

#### 3. GO Term Search
- Search by Gene Ontology term IDs (e.g., `GO:0003677`)
- Search by GO term names (e.g., "DNA binding", "transcription")
- Includes both experimental and computational annotations


## File Structure

```
plasmofp website/
├── plasmoFP_explorer_simple.py       # Streamlit web application
├── preprocess_enhanced_data.py       # Data preprocessing script
├── requirements.txt                  # Python dependencies
├── README.md                         # This file
├── optimized_gene_index.json         # Pre-processed gene index
├── search_index.json                 # Search optimization index
├── go_terms.json                     # Gene Ontology terms
├── protein_descriptions.json         # Protein description mappings
├── product_to_transcripts.json       # Product to gene mappings
├── product_search_index.json         # Product search index
├── transcript_to_species.json        # Gene to species mappings
├── go_id_to_genes.json               # GO ID to gene mappings
├── go_name_to_genes.json             # GO name to gene mappings
├── go_search_index.json              # GO term search index
├── annotated_fasta/                  # FASTA files directory
│   └── ... (19 species FASTA files)
└── with_PFP_predictions_complete_2/  # PlasmoFP data directory
    ├── PlasmoDB-68_Pfalciparum3D7_gene_dict_with_PFP.pkl
    ├── PlasmoDB-68_PvivaxSal1_gene_dict_with_PFP.pkl
    └── ... (17 more species files)
```



