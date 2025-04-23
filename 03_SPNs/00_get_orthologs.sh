#!/bin/bash
#module load python/3.8.x-anaconda
#conda activate /home2/gkonop/my_conda_env
#conda env remove --prefix /project/Neuroinformatics_Core/Konopka_lab/shared/For_Gozde/03_INTEGRATED_ALL/SPN/my_conda_env
#!/bin/bash

# Define input and output directories
input_file="/endosome/work/Neuroinformatics_Core/gkonop/03_INTEGRATE_ALL/SPN/human_pr_coding_genes.txt"
output_dir="/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/orthologs_output_json_2"
temp_dir="/home2/gkonop/projectShared/For_Gozde/03_INTEGRATED_ALL/SPN/temp_orthologs_2"

# Define taxIDs for the species of interest
declare -A species_taxids=(
    ["Homo sapiens"]=9606
    ["Pan troglodytes"]=9598
    ["Macaca mulatta"]=9544
    ["Callithrix jacchus"]=9483
    ["Mus musculus"]=10090
    ["Phyllostomus discolor"]=89673
    ["Mustela putorius furo"]=9669
)

# Create necessary directories (Uncomment if you need to create them)
#mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# Loop over each gene symbol in the input file
while IFS= read -r gene_id; do
    # Define the path to the .jsonl file for this gene
    jsonl_file="$output_dir/$gene_id.jsonl"

    # Check if the .jsonl file already exists
    if [ -f "$jsonl_file" ]; then
        echo "Gene $gene_id already processed. Skipping..." >> "$output_dir/already_processed.txt"
        continue  # Skip to the next gene
    else
        echo "Processing gene: $gene_id..."
    fi

    # Change to temp_dir before downloading
    cd "$temp_dir" || { echo "Failed to change directory to $temp_dir"; exit 1; }

    # Initialize an empty JSON array in the gene's .jsonl file if it doesn't exist yet
    echo "[]" > "$jsonl_file"  # Start with an empty JSON array for new genes

    # Loop over each species_taxid and download the ortholog data for each
    for species in "${!species_taxids[@]}"; do
        taxid="${species_taxids[$species]}"

        echo "Downloading orthologs for $gene_id from species $species (TaxID: $taxid)..."

        # Run the 'datasets' command to download the ortholog data for the current species
        /home2/gkonop/workdir/PROGRAMS/datasets download gene symbol "$gene_id" --ortholog "$taxid"

        # Check if the command completed successfully
        if [ $? -ne 0 ]; then
            echo "Error downloading data for $gene_id from species $species ($taxid). Skipping..." >> "$temp_dir/error_log.txt"
            continue
        fi

        # Check if the zip file exists (assuming it's downloaded into the current directory)
        zip_file=$(ls "$temp_dir"/*.zip)

        if [ ! -f "$zip_file" ]; then
            echo "Error: No zip file found for $gene_id from species $species ($taxid). Skipping..." >> "$temp_dir/error_log.txt"
            continue
        fi

        # Extract the .zip file to a temporary directory
        unzip -o "$zip_file" -d "$temp_dir"

        # Check if the data_report.jsonl file exists
        if [ ! -f "$temp_dir/ncbi_dataset/data/data_report.jsonl" ]; then
            echo "data_report.jsonl not found for $gene_id from species $species ($taxid). Skipping..." >> "$temp_dir/error_log.txt"
            continue
        fi

        # Append the content of data_report.jsonl to the gene's jsonl file
        # Read the contents of data_report.jsonl and append it as a single JSON object for the current species
        jq -c ".[] | {gene_id: \"$gene_id\", species: \"$species\", taxid: \"$taxid\", ortholog: .}" "$temp_dir/ncbi_dataset/data/data_report.jsonl" >> "$jsonl_file"

        # Clean up the extracted files to save space
        rm -rf "$temp_dir/ncbi_dataset"

        echo "Successfully processed gene $gene_id for species $species ($taxid)."
    done

done < "$input_file"

echo "Orthologs data download and extraction complete."
