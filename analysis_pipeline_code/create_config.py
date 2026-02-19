import yaml
from pathlib import Path

# 建立 config.yaml (Web views.py 和 MME_Web_Tool_CLI 共用)
def create_config(
    job_dir: Path,
    matching_logic: str,
    query_fasta_in_job: Path,
    human_fasta=None,
    human_db=None,
    IEDB_file=None,
    human_protein_detail=None,
    perfect_match_params=None,
    blastp_params=None,
):
    config = {
        "Matching_Logic": matching_logic,   # "Perfect_Match" or "Blastp"

        "perfect_match_paths": {
            "human_fasta_path": str(human_fasta),
            "query_fasta_path": str(query_fasta_in_job),
            "IEDB_file_path": str(IEDB_file),
            "human_protein_detail_path": str(human_protein_detail),
        },

        "perfect_match_parameters": perfect_match_params or {},

        "blastp_paths": {
            "human_db_path": str(human_db),
            "query_fasta_path": str(query_fasta_in_job),
            "IEDB_file_path": str(IEDB_file),
            "human_protein_detail_path": str(human_protein_detail),
        },

        "blastp_parameters": blastp_params or {},
    }

    output_path = job_dir / "config.yaml"

    with open(output_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(config, f, allow_unicode=True, sort_keys=False)

    print(f"📝 config.yaml 建立完成 → {output_path}")
