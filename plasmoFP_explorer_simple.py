#!/usr/bin/env python3
"""
PlasmoFP Explorer
"""

import streamlit as st
import json
import pandas as pd
import plotly.express as px
from collections import defaultdict
import uuid

# Page configuration
st.set_page_config(
    page_title="PlasmoFP Explorer",
    page_icon="DNA",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
    }
    
    .gene-header {
        background: linear-gradient(90deg, #1e3a8a 0%, #3b82f6 100%);
        color: white;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
    }
    
    .ontology-section {
        border: 1px solid #e5e7eb;
        border-radius: 0.5rem;
        padding: 1rem;
        margin-bottom: 1rem;
        background-color: #fafafa;
    }
    
    .prediction-score {
        background-color: #dbeafe;
        padding: 0.25rem 0.5rem;
        border-radius: 0.25rem;
        font-size: 0.875rem;
        font-weight: 600;
    }
    
    .go-term {
        background-color: #f3f4f6;
        padding: 0.25rem 0.5rem;
        border-radius: 0.25rem;
        margin: 0.125rem;
        display: inline-block;
        font-size: 0.875rem;
    }
    
    .stats-metric {
        text-align: center;
        background-color: #f9fafb;
        padding: 1rem;
        border-radius: 0.5rem;
        border: 1px solid #e5e7eb;
    }
</style>
""", unsafe_allow_html=True)

@st.cache_data
def load_cluster_mappings():
    """Load GO term cluster mappings from TSV files"""
    cluster_mappings = {}
    
    cluster_files = {
        'MF': 'MF_term_clustersAUG23.tsv',
        'BP': 'BP_term_clustersAUG23.tsv', 
        'CC': 'CC_term_clustersAUG23.tsv'
    }
    
    for aspect, filename in cluster_files.items():
        try:
            df = pd.read_csv(filename, sep='\t')
            cluster_mappings[aspect] = {}
            
            for _, row in df.iterrows():
                go_id = row['GO_ID']
                cluster_id = row['ClusterID']
                cluster_name = row['ClusterName']
                
                cluster_mappings[aspect][go_id] = {
                    'cluster_id': cluster_id,
                    'cluster_name': cluster_name
                }
        except FileNotFoundError:
            st.warning(f"Cluster file {filename} not found. Cluster analysis will be disabled for {aspect}.")
            cluster_mappings[aspect] = {}
        except Exception as e:
            st.error(f"Error loading cluster file {filename}: {e}")
            cluster_mappings[aspect] = {}
    
    return cluster_mappings

@st.cache_data
def load_optimized_data():
    """Load pre-processed optimized data files including enhanced search indices"""
    try:
        # Load gene index
        with open('optimized_gene_index.json', 'r') as f:
            gene_index = json.load(f)
        
        # Load search index
        with open('search_index.json', 'r') as f:
            search_index = json.load(f)
        
        # Load GO terms for fallback name resolution
        with open('go_terms.json', 'r') as f:
            go_terms = json.load(f)
        
        # Load cluster mappings
        cluster_mappings = load_cluster_mappings()
        
        # Load enhanced search indices
        with open('protein_descriptions.json', 'r') as f:
            protein_descriptions = json.load(f)
            
        with open('product_to_transcripts.json', 'r') as f:
            product_to_transcripts = json.load(f)
            
        with open('product_search_index.json', 'r') as f:
            product_search_index = json.load(f)
            
        with open('transcript_to_species.json', 'r') as f:
            transcript_to_species = json.load(f)
            
        with open('go_id_to_genes.json', 'r') as f:
            go_id_to_genes = json.load(f)
            
        with open('go_name_to_genes.json', 'r') as f:
            go_name_to_genes = json.load(f)
            
        with open('go_search_index.json', 'r') as f:
            go_search_index = json.load(f)
        
        enhanced_indices = {
            'protein_descriptions': protein_descriptions,
            'product_to_transcripts': product_to_transcripts,
            'product_search_index': product_search_index,
            'transcript_to_species': transcript_to_species,
            'go_id_to_genes': go_id_to_genes,
            'go_name_to_genes': go_name_to_genes,
            'go_search_index': go_search_index
        }
        
        return gene_index, search_index, go_terms, cluster_mappings, enhanced_indices
    
    except FileNotFoundError as e:
        st.error(f"Optimized data files not found: {e}")
        st.info("Please run 'python preprocess_enhanced_data.py' first to generate the enhanced search indices.")
        st.stop()
    except Exception as e:
        st.error(f"Error loading data: {e}")
        st.stop()

def search_genes(query, gene_index, search_index, limit=50):
    """Fast gene search using pre-built index. Supports comma-separated multi-gene search."""
    if not query:
        return []
    
    # Check if this is a multi-gene search (contains commas)
    if ',' in query:
        return search_multiple_genes(query, gene_index)
    
    # Single gene search logic
    query_upper = query.upper()
    candidate_genes = set()
    
    # Use search index for fast prefix matching
    for prefix, genes in search_index.items():
        if query_upper in prefix:
            candidate_genes.update(genes)
    
    # If no candidates from index, fall back to direct search
    if not candidate_genes:
        for gene_id in gene_index.keys():
            if query_upper in gene_id.upper():
                candidate_genes.add(gene_id)
    
    # Filter and limit results
    results = []
    for gene_id in candidate_genes:
        if gene_id in gene_index and query_upper in gene_id.upper():
            results.append((gene_id, gene_index[gene_id]))
            if len(results) >= limit:
                break
    
    # Sort results to prioritize exact matches
    results.sort(key=lambda x: (
        not x[0].upper().startswith(query_upper),  # Exact prefix matches first
        len(x[0]),  # Shorter gene IDs first
        x[0]  # Alphabetical order
    ))
    
    return results

def search_multiple_genes(query, gene_index):
    """Search for multiple genes separated by commas. Returns exact matches only."""
    # Split by comma and clean up each gene ID
    gene_ids = [gene_id.strip() for gene_id in query.split(',')]
    
    results = []
    for gene_id in gene_ids:
        if not gene_id:  # Skip empty strings
            continue
            
        # Try exact match first
        if gene_id in gene_index:
            results.append((gene_id, gene_index[gene_id]))
        else:
            # Try case-insensitive exact match
            found = False
            for existing_gene_id in gene_index.keys():
                if existing_gene_id.upper() == gene_id.upper():
                    results.append((existing_gene_id, gene_index[existing_gene_id]))
                    found = True
                    break
            
            # If still not found, try partial matching for this specific gene
            if not found:
                gene_id_upper = gene_id.upper()
                for existing_gene_id in gene_index.keys():
                    if gene_id_upper in existing_gene_id.upper():
                        results.append((existing_gene_id, gene_index[existing_gene_id]))
                        found = True
                        break
    
    return results

def search_genes_by_product(query, enhanced_indices, gene_index, limit=1000):
    """Search genes by protein description with fuzzy matching"""
    if not query:
        return []
    
    query_lower = query.lower()
    words = query_lower.split()
    
    # Find matching products using the search index
    matching_products = set()
    
    for word in words:
        # Find products that contain this word (or prefixes)
        for search_term, products in enhanced_indices['product_search_index'].items():
            if word in search_term or search_term in word:
                matching_products.update(products)
    
    # If no fuzzy matches, try direct substring search
    if not matching_products:
        for product in enhanced_indices['product_to_transcripts'].keys():
            if query_lower in product.lower():
                matching_products.add(product)
    
    # Get genes for matching products and group by species
    results_by_species = defaultdict(list)
    total_genes = 0
    
    for product in matching_products:
        genes_with_product = enhanced_indices['product_to_transcripts'][product]
        
        for gene_id in genes_with_product:
            if total_genes >= limit:
                break
                
            # Check if gene exists in our gene index
            if gene_id in gene_index:
                gene_data = gene_index[gene_id]
                species = enhanced_indices['transcript_to_species'].get(gene_id, gene_data.get('species', 'Unknown'))
                
                results_by_species[species].append({
                    'gene_id': gene_id,
                    'gene_data': gene_data,
                    'product': product
                })
                total_genes += 1
        
        if total_genes >= limit:
            break
    
    return results_by_species

def search_genes_by_go_term(query, enhanced_indices, gene_index, go_terms, limit=1000):
    """Search genes by GO term ID or name"""
    if not query:
        return []
    
    query_lower = query.lower().strip()
    matching_genes = set()
    
    # Check if query is a GO ID (starts with GO:)
    if query_lower.startswith('go:'):
        go_id = query_lower.upper()
        if go_id in enhanced_indices['go_id_to_genes']:
            for gene_id, species, source in enhanced_indices['go_id_to_genes'][go_id]:
                matching_genes.add((gene_id, species, source, go_id, go_terms.get(go_id, go_id)))
    else:
        # Search by GO term name
        words = query_lower.split()
        matching_go_names = set()
        
        # Find matching GO names using fuzzy search
        for word in words:
            for search_term, go_names in enhanced_indices['go_search_index'].items():
                if word in search_term or search_term in word:
                    matching_go_names.update(go_names)
        
        # Also try direct substring search
        for go_name in enhanced_indices['go_name_to_genes'].keys():
            if query_lower in go_name:
                matching_go_names.add(go_name)
        
        # Get genes for matching GO names
        for go_name in matching_go_names:
            if go_name in enhanced_indices['go_name_to_genes']:
                for gene_id, species, source in enhanced_indices['go_name_to_genes'][go_name]:
                    matching_genes.add((gene_id, species, source, 'N/A', go_name))
    
    # Group results by species and limit
    results_by_species = defaultdict(list)
    total_genes = 0
    
    for gene_id, species, source, go_id, go_name in matching_genes:
        if total_genes >= limit:
            break
            
        if gene_id in gene_index:
            gene_data = gene_index[gene_id]
            results_by_species[species].append({
                'gene_id': gene_id,
                'gene_data': gene_data,
                'go_id': go_id,
                'go_name': go_name,
                'source': source
            })
            total_genes += 1
    
    return results_by_species

def create_annotation_source_pie_chart(results_by_species):
    """Create pie chart showing distribution of annotation sources"""
    source_counts = defaultdict(int)
    
    for species_results in results_by_species.values():
        for result in species_results:
            source = result.get('source', 'Unknown')
            source_counts[source] += 1
    
    if not source_counts:
        return None
    
    # Create pie chart
    fig = px.pie(
        values=list(source_counts.values()),
        names=list(source_counts.keys()),
        title="Annotation Source Distribution"
    )
    
    fig.update_traces(
        textposition='inside',
        textinfo='percent+label'
    )
    
    fig.update_layout(
        showlegend=True,
        height=400,
        margin=dict(t=50, b=50, l=50, r=50)
    )
    
    return fig

def create_species_distribution_pie_chart(results_by_species):
    """Create pie chart showing distribution across species"""
    species_counts = {species: len(results) for species, results in results_by_species.items()}
    
    if not species_counts:
        return None
    
    # Create pie chart
    fig = px.pie(
        values=list(species_counts.values()),
        names=list(species_counts.keys()),
        title="Species Distribution"
    )
    
    fig.update_traces(
        textposition='inside',
        textinfo='percent+label'
    )
    
    fig.update_layout(
        showlegend=True,
        height=400,
        margin=dict(t=50, b=50, l=50, r=50)
    )
    
    return fig

def create_cluster_distribution_chart(predictions, aspect_code, cluster_mappings):
    """Create pie chart showing cluster distribution for predictions"""
    if not predictions or aspect_code not in cluster_mappings:
        return None
    
    # Count terms by cluster
    cluster_counts = defaultdict(int)
    cluster_names = {}
    
    for pred in predictions:
        go_id = pred['id']
        if go_id in cluster_mappings[aspect_code]:
            cluster_info = cluster_mappings[aspect_code][go_id]
            cluster_id = cluster_info['cluster_id']
            cluster_name = cluster_info['cluster_name']
            
            cluster_counts[cluster_id] += 1
            cluster_names[cluster_id] = cluster_name
    
    if not cluster_counts:
        return None
    
    # Prepare data for pie chart
    chart_data = []
    for cluster_id, count in cluster_counts.items():
        chart_data.append({
            'Cluster': f"{cluster_names[cluster_id]}",
            'Count': count,
            'Cluster_ID': cluster_id
        })
    
    df = pd.DataFrame(chart_data)
    
    # Create pie chart
    fig = px.pie(
        df, 
        values='Count', 
        names='Cluster',
        title=f"Functional Cluster Distribution",
        hover_data=['Cluster_ID']
    )
    
    # Update layout for better appearance
    fig.update_traces(
        textposition='inside', 
        textinfo='percent+label',
        hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<br>Cluster ID: %{customdata[0]}<extra></extra>'
    )
    
    fig.update_layout(
        showlegend=True,
        height=400,
        margin=dict(t=50, b=50, l=50, r=50)
    )
    
    return fig

def display_ontology_data(ontology_name, icon, plasmofp_data, original_data, selected_fdr, go_terms=None, cluster_mappings=None, aspect_code=None, gene_id=None, gene_data=None):
    """Display data for a specific ontology (MF, BP, CC)"""
    st.markdown(f"### {icon} {ontology_name}")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### PlasmoFP Predictions")
        # Convert selected_fdr to string to match data format
        fdr_key = str(selected_fdr)
        if plasmofp_data and fdr_key in plasmofp_data:
            predictions = plasmofp_data[fdr_key]
            if predictions:
                st.success(f"Found {len(predictions)} predictions at {selected_fdr:.2f} eFDR")
                
                # Create DataFrame for predictions with improved GO term name resolution and cluster info
                pred_data = []
                for pred in predictions:
                    go_id = pred['id']
                    go_name = pred['name']
                    score = pred['score']
                    
                    # If the name is just the GO ID, try to get the actual name from GO terms
                    if go_name == go_id and go_terms and go_id in go_terms:
                        go_name = go_terms[go_id]
                    elif go_name == go_id:
                        # If still no name found, add a note about missing term
                        go_name = f"{go_id} (term not in current GO release)"
                    
                    # Get cluster information
                    cluster_info = "(not present in clustering file)"
                    if cluster_mappings and aspect_code in cluster_mappings and go_id in cluster_mappings[aspect_code]:
                        cluster_data = cluster_mappings[aspect_code][go_id]
                        cluster_info = f"{cluster_data['cluster_id']}: {cluster_data['cluster_name']}"
                    
                    pred_data.append({
                        'GO ID': go_id,
                        'GO Term': go_name,
                        'Score': round(score, 4),
                        'Cluster ID': cluster_info
                    })
                
                pred_df = pd.DataFrame(pred_data)
                pred_df = pred_df.sort_values('Score', ascending=False)
                
                st.dataframe(pred_df, use_container_width=True, hide_index=True)
            else:
                st.info(f"No predictions at {selected_fdr:.2f} eFDR")
        else:
            st.info("No PlasmoFP predictions available")
    
    with col2:
        st.markdown("#### Original Annotations")
        if original_data:
            # Deduplicate original annotations based on GO ID
            seen_go_ids = set()
            unique_annotations = []
            for annotation in original_data:
                go_id = annotation['id']
                if go_id not in seen_go_ids:
                    unique_annotations.append(annotation)
                    seen_go_ids.add(go_id)
            
            st.success(f"Found {len(unique_annotations)} original annotations")
            
            # Create DataFrame for original annotations with improved GO term name resolution and cluster info
            orig_data_processed = []
            for annotation in unique_annotations:
                go_id = annotation['id']
                go_name = annotation['name']
                
                # If the name is just the GO ID, try to get the actual name from GO terms
                if go_name == go_id and go_terms and go_id in go_terms:
                    go_name = go_terms[go_id]
                elif go_name == go_id:
                    # If still no name found, add a note about missing term
                    go_name = f"{go_id} (term not in current GO release)"
                
                # Get cluster information
                cluster_info = "(not present in clustering file)"
                if cluster_mappings and aspect_code in cluster_mappings and go_id in cluster_mappings[aspect_code]:
                    cluster_data = cluster_mappings[aspect_code][go_id]
                    cluster_info = f"{cluster_data['cluster_id']}: {cluster_data['cluster_name']}"
                
                orig_data_processed.append({
                    'GO ID': go_id,
                    'GO Term': go_name,
                    'Cluster ID': cluster_info
                })
            
            orig_df = pd.DataFrame(orig_data_processed)
            st.dataframe(orig_df, use_container_width=True, hide_index=True)
        else:
            st.info("No original annotations available")
    
    # Add cluster distribution chart below the tables
    if plasmofp_data and fdr_key in plasmofp_data and plasmofp_data[fdr_key]:
        st.markdown("#### Functional Cluster Distribution")
        cluster_chart = create_cluster_distribution_chart(plasmofp_data[fdr_key], aspect_code, cluster_mappings)
        
        if cluster_chart:
            st.plotly_chart(cluster_chart, use_container_width=True, key=f"cluster_{gene_id.replace('.', '_').replace('-', '_')}_{aspect_code}_{str(selected_fdr).replace('.', '')}")
        else:
            st.info("No clustered terms available for visualization")

def display_gene_info(gene_id, gene_data, go_terms=None, cluster_mappings=None):
    """Display gene information in the new simplified format"""
    
    # Gene header
    st.markdown(f"""
    <div class="gene-header">
        <h2>Gene: {gene_id}</h2>
        <p><strong>Species:</strong> {gene_data['species']}</p>
    </div>
    """, unsafe_allow_html=True)
    
    # eFDR threshold selector
    st.markdown("#### eFDR Threshold Selection")
    
    # Get all available eFDR thresholds
    all_fdrs = set()
    for aspect_data in gene_data.get('pfp_predictions', {}).values():
        all_fdrs.update(aspect_data.keys())
    
    if all_fdrs:
        # Convert to floats and sort
        fdr_options = sorted([float(fdr) for fdr in all_fdrs])
        
        # Use a simple, stable key based on gene ID only
        selected_fdr = st.selectbox(
            "Select eFDR threshold:",
            fdr_options,
            index=1 if len(fdr_options) > 1 else 0,
            format_func=lambda x: f"{x:.2f} ({x*100:.0f}%)",
            key=f"efdr_{gene_id.replace('.', '_').replace('-', '_')}"
        )
    else:
        selected_fdr = 0.05  # Default fallback
        st.info("No PlasmoFP predictions available - showing original annotations only")
    
    st.markdown("---")
    
    # Display data by ontology
    ontologies = [
        ("Molecular Function", "MF", "MF"),
        ("Biological Process", "BP", "BP"), 
        ("Cellular Component", "CC", "CC")
    ]
    
    for ontology_name, icon, aspect_code in ontologies:
        # Get data for this ontology
        plasmofp_data = gene_data.get('pfp_predictions', {}).get(aspect_code, {})
        original_data = gene_data.get('original_annotations', {}).get(aspect_code, [])
        
        # Always show all ontology sections for consistency
        with st.container():
            st.markdown(f'<div class="ontology-section">', unsafe_allow_html=True)
            display_ontology_data(ontology_name, icon, plasmofp_data, original_data, selected_fdr, go_terms, cluster_mappings, aspect_code, gene_id, gene_data)
            st.markdown('</div>', unsafe_allow_html=True)

def main():
    """Main application function"""
    
    st.title("PlasmoFP Explorer")
    # st.markdown("**Explore Protein Function Predictions for Plasmodium Species**")
    
    # About section - always visible at top
    st.markdown("""
    **PlasmoFP Explorer** is a Streamlit web application for exploring protein function prediction 
    for *Plasmodium* species. Predictions are generated by PlasmoFP, a structure-informed deep learning model developed specifically for *Plasmodium* species. For more infromation on PlasmoFP, please refer to the manuscript: PlasmoFP: leveraging deep learning to predict protein function of uncharacterized proteins across the malaria parasite genus 
    
    **Features:**
    - **Gene ID Search**: Look up genes by their identifiers
    - **Protein Description Search**: Find genes by their product descriptions
    - **GO Term Search**: Search by Gene Ontology terms or Gene Ontology term descriptors 
    - **Multi-species Coverage**: 19 *Plasmodium* species 
    - **Functional Clustering**: Terms organized by functional clusters 

    Developed by Jane Carlton's lab at Johns Hopkins University Malaria Research Institute
    """)
    
    st.markdown("---")
    
    # Load data
    with st.spinner("Loading data..."):
        gene_index, search_index, go_terms, cluster_mappings, enhanced_indices = load_optimized_data()
    
    # Sidebar stats
    st.sidebar.title("Database Stats")
    
    total_genes = len(gene_index)
    genes_with_plasmofp = sum(1 for data in gene_index.values() if data.get('pfp_predictions'))
    genes_with_original = sum(1 for data in gene_index.values() if data.get('original_annotations'))
    
    # Species counts
    species_counts = defaultdict(int)
    for data in gene_index.values():
        species_counts[data['species']] += 1
    
    st.sidebar.metric("Total Genes", f"{total_genes:,}")
    st.sidebar.metric("Genes with PlasmoFP", f"{genes_with_plasmofp:,}")
    st.sidebar.metric("Genes with Annotations", f"{genes_with_original:,}")
    st.sidebar.metric("Species", len(species_counts))
    
    # Main search interface
    st.markdown("## Search Interface")
    
    # Search type selector
    search_type = st.radio(
        "Search by:",
        ["Gene ID", "Protein Description", "GO Term"],
        horizontal=True,
        key="search_type"
    )
    
    # Dynamic placeholders based on search type
    placeholders = {
        "Gene ID": "e.g., PF3D7_1206100, PVX_000005, PmUG01_... or comma-separated for multiple",
        "Protein Description": "e.g., unknown function, hypothetical protein, DNA helicase, protein kinase",
        "GO Term": "e.g., GO:0003677, DNA binding, transcription, kinase activity"
    }
    
    col1, col2 = st.columns([4, 1])
    
    with col1:
        query = st.text_input(
            f"Enter {search_type.lower()}:",
            placeholder=placeholders[search_type],
            key="main_search"
        )
    
    with col2:
        max_results = st.selectbox("Max results:", [25, 50, 100, 250, 500], index=2)
    
    # Search and display results
    if query:
        with st.spinner("Searching..."):
            if search_type == "Gene ID":
                # Original gene ID search
                results = search_genes(query, gene_index, search_index, limit=max_results)
                
                if results:
                    # Check if this was a multi-gene search
                    is_multi_search = ',' in query
                    if is_multi_search:
                        requested_genes = [g.strip() for g in query.split(',') if g.strip()]
                        st.success(f"Found {len(results)} of {len(requested_genes)} requested genes")
                    else:
                        st.success(f"Found {len(results)} matching genes")
                    
                    # Show results in expandable sections
                    for i, (gene_id, gene_data) in enumerate(results):
                        with st.expander(
                            f"{gene_id} ({gene_data['species']})", 
                            expanded=(i == 0 and len(results) <= 3)
                        ):
                                display_gene_info(gene_id, gene_data, go_terms, cluster_mappings)
                else:
                    st.warning("No genes found matching your search criteria.")
                    st.info("Try a different search term or check the spelling.")
                    
            elif search_type == "Protein Description":
                # Protein description search
                results_by_species = search_genes_by_product(query, enhanced_indices, gene_index, limit=max_results)
                
                if results_by_species:
                    total_results = sum(len(results) for results in results_by_species.values())
                    if total_results >= max_results:
                        st.success(f"Showing {total_results} genes with matching protein descriptions across {len(results_by_species)} species")
                        st.info(f"Results limited to {max_results}. Increase 'Max results' to see more.")
                    else:
                        st.success(f"Found {total_results} genes with matching protein descriptions across {len(results_by_species)} species")
                    
                    # Display results grouped by species
                    for species, results in results_by_species.items():
                        st.subheader(f"Species: {species} ({len(results)} genes)")
                        
                        for i, result in enumerate(results):
                            gene_id = result['gene_id']
                            gene_data = result['gene_data']
                            product = result['product']
                            
                            # Truncate long product descriptions for display
                            display_product = product[:100] + "..." if len(product) > 100 else product
                            
                            with st.expander(
                                f"{gene_id} - {display_product}", 
                                expanded=False
                            ):
                                st.markdown(f"**Full Product Description:** {product}")
                                display_gene_info(gene_id, gene_data, go_terms, cluster_mappings)
                else:
                    st.warning("No genes found with matching protein descriptions.")
                    st.info("Try different keywords or check the spelling.")
                    
            elif search_type == "GO Term":
                # GO term search
                results_by_species = search_genes_by_go_term(query, enhanced_indices, gene_index, go_terms, limit=max_results)
                
                if results_by_species:
                    total_results = sum(len(results) for results in results_by_species.values())
                    if total_results >= max_results:
                        st.success(f"Showing {total_results} genes with matching GO terms across {len(results_by_species)} species")
                        st.info(f"Results limited to {max_results}. Increase 'Max results' to see more.")
                    else:
                        st.success(f"Found {total_results} genes with matching GO terms across {len(results_by_species)} species")
                    
                    # Display results grouped by species
                    for species, results in results_by_species.items():
                        st.subheader(f"Species: {species} ({len(results)} genes)")
                        
                        for i, result in enumerate(results):
                            gene_id = result['gene_id']
                            gene_data = result['gene_data']
                            go_name = result['go_name']
                            source = result['source']
                            
                            with st.expander(
                                f"{gene_id} - {go_name} ({source})", 
                                expanded=False
                            ):
                                display_gene_info(gene_id, gene_data, go_terms, cluster_mappings)
                    
                    # Add summary pie charts for GO term results
                    st.markdown("---")
                    st.markdown("### Summary Visualizations")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        source_chart = create_annotation_source_pie_chart(results_by_species)
                        if source_chart:
                            st.plotly_chart(source_chart, use_container_width=True, key=f"go_source_chart_{hash(query + 'source')}")
                    
                    with col2:
                        species_chart = create_species_distribution_pie_chart(results_by_species)
                        if species_chart:
                            st.plotly_chart(species_chart, use_container_width=True, key=f"go_species_chart_{hash(query + 'species')}")
                            
                else:
                    st.warning("No genes found with matching GO terms.")
                    st.info("Try different GO terms, IDs, or check the spelling.")
    

if __name__ == "__main__":
    main()
