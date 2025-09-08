#!/usr/bin/env python3
"""
Enhanced Data Preprocessing for PlasmoFP Explorer
Parses FASTA files and creates optimized search indices for:
1. Gene ID search (existing)
2. Protein description search (new)
3. GO term search (new)
"""

import json
import pickle
import pandas as pd
from collections import defaultdict
import re
import os
from pathlib import Path

def parse_fasta_header(header_line):
    """Parse FASTA header to extract gene information"""
    # Remove the '>' and split by ' | '
    parts = header_line[1:].split(' | ')
    
    gene_info = {}
    for part in parts:
        if '=' in part:
            key, value = part.split('=', 1)
            gene_info[key] = value
        else:
            # First part is the gene ID
            gene_info['gene_id'] = part
    
    return gene_info

def load_fasta_files(fasta_dir):
    """Load and parse all FASTA files to extract gene product mappings"""
    print("Loading FASTA files...")
    
    gene_to_product = {}
    product_to_genes = defaultdict(list)
    gene_to_species = {}
    
    fasta_files = list(Path(fasta_dir).glob("*.fasta"))
    
    for fasta_file in fasta_files:
        print(f"Processing {fasta_file.name}...")
        
        with open(fasta_file, 'r') as f:
            for line_num, line in enumerate(f):
                if line.startswith('>'):
                    try:
                        gene_info = parse_fasta_header(line.strip())
                        
                        if 'gene_id' in gene_info and 'gene_product' in gene_info:
                            gene_id = gene_info['gene_id']  # Use the actual gene ID from header
                            product = gene_info['gene_product']
                            organism = gene_info.get('organism', 'Unknown')
                            
                            # Extract species name from organism
                            species = organism.replace('_', ' ').replace('Plasmodium ', 'P')
                            
                            gene_to_product[gene_id] = product
                            product_to_genes[product].append(gene_id)
                            gene_to_species[gene_id] = species
                            
                    except Exception as e:
                        print(f"Error parsing line {line_num} in {fasta_file.name}: {e}")
                        continue
    
    print(f"Loaded {len(gene_to_product)} gene-product mappings")
    return gene_to_product, product_to_genes, gene_to_species

def create_product_search_index(product_to_genes):
    """Create search index for fuzzy protein description matching"""
    print("Creating product search index...")
    
    search_index = defaultdict(set)
    
    for product in product_to_genes.keys():
        # Split product into words and create searchable terms
        words = re.findall(r'\b\w+\b', product.lower())
        
        for word in words:
            if len(word) >= 3:  # Only index words with 3+ characters
                # Add the word itself
                search_index[word].add(product)
                
                # Add partial matches (prefixes)
                for i in range(3, len(word) + 1):
                    prefix = word[:i]
                    search_index[prefix].add(product)
    
    # Convert sets to lists for JSON serialization
    search_index = {k: list(v) for k, v in search_index.items()}
    
    print(f"Created search index with {len(search_index)} terms")
    return search_index

def extract_go_mappings(gene_index):
    """Extract GO term to gene mappings from existing gene data"""
    print("Extracting GO term mappings...")
    
    go_id_to_genes = defaultdict(set)
    go_name_to_genes = defaultdict(set)
    
    for gene_id, gene_data in gene_index.items():
        species = gene_data.get('species', 'Unknown')
        
        # Process PlasmoFP predictions
        pfp_predictions = gene_data.get('pfp_predictions', {})
        for aspect in ['MF', 'BP', 'CC']:
            if aspect in pfp_predictions:
                for fdr_threshold, predictions in pfp_predictions[aspect].items():
                    for pred in predictions:
                        go_id = pred['id']
                        go_name = pred['name']
                        
                        go_id_to_genes[go_id].add((gene_id, species, 'PlasmoFP'))
                        
                        # Also index by GO name for search
                        if go_name and go_name != go_id:
                            go_name_to_genes[go_name.lower()].add((gene_id, species, 'PlasmoFP'))
        
        # Process original annotations
        original_annotations = gene_data.get('original_annotations', {})
        for aspect in ['MF', 'BP', 'CC']:
            if aspect in original_annotations:
                for annotation in original_annotations[aspect]:
                    go_id = annotation['id']
                    go_name = annotation['name']
                    
                    go_id_to_genes[go_id].add((gene_id, species, 'Original'))
                    
                    if go_name and go_name != go_id:
                        go_name_to_genes[go_name.lower()].add((gene_id, species, 'Original'))
    
    # Convert sets to lists for JSON serialization
    go_id_to_genes = {k: list(v) for k, v in go_id_to_genes.items()}
    go_name_to_genes = {k: list(v) for k, v in go_name_to_genes.items()}
    
    print(f"Extracted mappings for {len(go_id_to_genes)} GO IDs and {len(go_name_to_genes)} GO names")
    return go_id_to_genes, go_name_to_genes

def create_go_search_index(go_name_to_genes):
    """Create fuzzy search index for GO term names"""
    print("Creating GO term search index...")
    
    search_index = defaultdict(set)
    
    for go_name in go_name_to_genes.keys():
        words = re.findall(r'\b\w+\b', go_name.lower())
        
        for word in words:
            if len(word) >= 3:
                search_index[word].add(go_name)
                
                # Add partial matches
                for i in range(3, len(word) + 1):
                    prefix = word[:i]
                    search_index[prefix].add(go_name)
    
    # Convert sets to lists
    search_index = {k: list(v) for k, v in search_index.items()}
    
    print(f"Created GO search index with {len(search_index)} terms")
    return search_index

def main():
    """Main preprocessing function"""
    print("Starting enhanced data preprocessing...")
    
    # Load existing optimized data
    print("Loading existing gene index...")
    with open('optimized_gene_index.json', 'r') as f:
        gene_index = json.load(f)
    
    # Load FASTA files for protein descriptions
    gene_to_product, product_to_genes, gene_to_species = load_fasta_files('annotated_fasta')
    
    # Create protein description search index
    product_search_index = create_product_search_index(product_to_genes)
    
    # Extract GO term mappings
    go_id_to_genes, go_name_to_genes = extract_go_mappings(gene_index)
    
    # Create GO term search index
    go_search_index = create_go_search_index(go_name_to_genes)
    
    # Save all new indices
    print("Saving enhanced search indices...")
    
    with open('protein_descriptions.json', 'w') as f:
        json.dump(gene_to_product, f)
    
    with open('product_to_transcripts.json', 'w') as f:
        json.dump(dict(product_to_genes), f)
    
    with open('product_search_index.json', 'w') as f:
        json.dump(product_search_index, f)
    
    with open('transcript_to_species.json', 'w') as f:
        json.dump(gene_to_species, f)
    
    with open('go_id_to_genes.json', 'w') as f:
        json.dump(go_id_to_genes, f)
    
    with open('go_name_to_genes.json', 'w') as f:
        json.dump(go_name_to_genes, f)
    
    with open('go_search_index.json', 'w') as f:
        json.dump(go_search_index, f)
    
    print("Enhanced preprocessing complete!")
    print(f"- Protein descriptions: {len(gene_to_product)}")
    print(f"- Unique products: {len(product_to_genes)}")
    print(f"- GO ID mappings: {len(go_id_to_genes)}")
    print(f"- GO name mappings: {len(go_name_to_genes)}")

if __name__ == "__main__":
    main()
