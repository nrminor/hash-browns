process LABKEY_UPLOAD_BLAST {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(blast_results), path(mapped_reads), path(virus_contigs)

    output:
    path "${sample_id}.blast.tkn"
    path "${sample_id}.csv"

    script:
    // experiment = config["global"]["experiment"],
    // blast_db_name = config["global"]["blast_db_name"]
    // snakemake_run_id = snakemake_run_id
    // labkey_server = config["global"].get("labkey_server" "")
    // project_name = config["global"].get("labkey_project_name" "")
    // api_key = config["global"].get("labkey_api_key" "")
    // stat_db_version = config["global"]["stat_dbss"]
    // out_dir = config["global"]["out_dir"]
    """
    labkey_upload.py
    """
}

process LABKEY_UPLOAD_FASTA {
    input:
    tuple val(sample_id), path(contigs)

    output:
    path "${sample_id}.fasta.tkn"

    script:
    // experiment = config["global"]["experiment"]
    // labkey_server = config["global"].get("labkey_server" "")
    // project_name = config["global"].get("labkey_project_name" "")        
    // api_key = config["global"].get("labkey_api_key" "")
    // snakemake_run_id = snakemake_run_id
    // singleline_fasta_path = f"{local_output_dir}/19_labkey/{{sample}}.fasta"
    // TODO: This script is hardcoded only to work with snakemake
    """
    labkey_upload_fasta.py
    """
}

process CREATE_TARBALL {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(virus_contigs), path(unclassified_contigs), path(mapped_reads), path(mapped_reads_index)

    output:
    path "${sample_id}.tar.zst"

    script:
    """
    # Initialize an empty list to hold the valid files
    valid_files=""

    # Check each file for existence and non-zero size
    for file in {input.virus_contigs} {input.unclassified_contigs} {input.mapped_reads} {input.mapped_reads_index}; do
        if [ -s "\$file" ]; then
            valid_files="\$valid_files \$file"
        else
            echo "Skipping missing or empty file: \$file"
        fi
    done

    # Only create the tarball if there are valid files
    if [ -n "\$valid_files" ]; then
        tar -I zstd -cvf {output.tarball} \$valid_files
    else
        echo "No valid files to include in the tarball."
    fi
    """
}

process UPLOAD_FILES_TO_LABKEY {
    input:
    path tarball

    script:
    // labkey_server = config["global"].get("labkey_server", ""),
    // project_name = config["global"].get("labkey_project_name", ""),   
    // webdav_url = config["global"].get("webdav_url", ""),
    // username = config["global"].get("labkey_username", ""),
    // password = config["global"].get("labkey_password", ""),
    // experiment = config["global"]["experiment"],
    // snakemake_run_id = snakemake_run_id,
    // final_output_path = config["global"]["out_dir"]

    // TODO: This script is hardcoded only to work with snakemake
    """
    upload_files_to_labkey.py
    """
}
