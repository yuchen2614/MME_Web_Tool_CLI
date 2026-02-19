from pathlib import Path
import uuid
import os
import yaml
import shutil

from MME_Web_Tool_CLI.analysis_pipeline_code.create_config import create_config
from MME_Web_Tool_CLI.analysis_pipeline_code.perfect_match_analysis_pipeline import perfect_match_pipeline_main
from MME_Web_Tool_CLI.analysis_pipeline_code.blastp_analysis_pipeline import blastp_pipeline_main


# command line 模式分析流程啟動程式
if __name__ == "__main__":

    # === Step 1 : 取得專案根目錄 ===
    project_root = Path(__file__).resolve().parent.parent
    jobs_root = project_root / "jobs"

    # === Step 2 : 生成不重複的JobID資料夾 (while 無限生成直到不重複) ===
    while True:
        job_id = uuid.uuid4().hex[:8]
        job_dir = jobs_root / job_id
        if not job_dir.exists():
            break   # ✅ 找到不重複 ID
    print(f"✅ 自動產生唯一 job_id: {job_id}")
    job_dir.mkdir(exist_ok=False)


    # === 正式分析流程 ===
    # === Step 3 : 讀取 CLI_config.yaml 的參數 ===
    with open(project_root / "analysis_pipeline_code" / "CLI_config.yaml", "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    logic = config.get("Matching_Logic")

    # 根據不同邏輯讀取不同參數
    if logic == "Perfect_Match":
        human_fasta_path = Path(config["perfect_match_paths"]["human_fasta_path"]).resolve()
        query_fasta_path = Path(config["perfect_match_paths"]["query_fasta_path"]).resolve()
        human_protein_detail_path = Path(config["perfect_match_paths"]["human_protein_detail_path"]).resolve()
        IEDB_file_path = Path(config["perfect_match_paths"]["IEDB_file_path"]).resolve()

        kmer_perfect_match = config["perfect_match_parameters"]["kmer_perfect_match"]

    elif logic == "Blastp":
        human_db_path = Path(config["blastp_paths"]["human_db_path"]).resolve()
        query_fasta_path = Path(config["blastp_paths"]["query_fasta_path"]).resolve()
        human_protein_detail_path = Path(config["blastp_paths"]["human_protein_detail_path"]).resolve()
        IEDB_file_path = Path(config["blastp_paths"]["IEDB_file_path"]).resolve()

        evalue = config["blastp_parameters"]["evalue"]
        matrix = config["blastp_parameters"]["matrix"]
        gap = config["blastp_parameters"]["gap"]
        kmer_blastp = config["blastp_parameters"]["kmer_blastp"]
        min_length = config["blastp_parameters"]["min_length"]
        max_length = config["blastp_parameters"]["max_length"]
        identity = config["blastp_parameters"]["identity"]
        positive = config["blastp_parameters"]["positive"]
        gap_num = config["blastp_parameters"]["gap_num"]

    # 修改 query fasta 的路徑，改為 job_id 資料夾，並複製一份 fasta 到 job_id 資料夾
    original_filename = query_fasta_path.name
    job_query_fasta_path = job_dir / original_filename
    shutil.copy(query_fasta_path, job_query_fasta_path)
    print(f"📁 Query FASTA 已複製到: {job_query_fasta_path}")

    # === Step 4 : 在 job_id 資料夾建立 config.yaml ===
    create_config(
        job_dir=job_dir,
        matching_logic=logic,

        query_fasta_in_job=job_query_fasta_path,
        human_fasta=human_fasta_path if logic == "Perfect_Match" else None,
        human_db=human_db_path if logic == "Blastp" else None,
        IEDB_file=IEDB_file_path,
        human_protein_detail=human_protein_detail_path,

        perfect_match_params={
            "kmer_perfect_match": kmer_perfect_match
        } if logic == "Perfect_Match" else None,

        blastp_params={
            "evalue": evalue,
            "matrix": matrix,
            "gap": gap,
            "kmer_blastp": kmer_blastp,
            "min_length": min_length,
            "max_length": max_length,
            "identity": identity,
            "positive": positive,
            "gap_num": gap_num,
        } if logic == "Blastp" else None
    )

    # === Step 5 : 執行正式分析 pipeline_main ===
    try: 
        if logic not in {"Perfect_Match", "Blastp"}:
            raise ValueError(f"Unknown Matching_Logic: {logic}")
        elif logic == "Perfect_Match":
            perfect_match_pipeline_main(job_id)
        elif logic == "Blastp":
            blastp_pipeline_main(job_id)
    except Exception as e:
        print("❌ Pipeline Main 發生錯誤：\n")