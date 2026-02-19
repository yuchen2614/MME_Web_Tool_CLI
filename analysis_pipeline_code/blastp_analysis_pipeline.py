import re
import sys
import yaml
import traceback
import subprocess
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from MME_Web_Tool_CLI.analysis_pipeline_code.shared_analysis_steps import (
    IEDB_analyze_matches,
    calc_reference_table,
    calc_epitope_table,
    calc_query_table,
)

#切片在Def加上座標，保留序列原始長度
def sliding_windows_fasta(fasta_text, window_size, group_size=100):
    """
    對多條 FASTA 序列進行 sliding window 切片，每 group_size 條分組。

    回傳:
        fasta_groups: List[str]，分組好的 FASTA 序列
        original_lengths: Dict[str, int]，記錄每條原始序列長度
    """
    fasta_groups = []
    current_group = []
    original_lengths = {}

    lines = fasta_text.strip().splitlines()
    title = None
    sequence = ""

    def process_sequence(title, sequence):
        result = []
        if not title or not sequence:
            return result
        original_lengths[title] = len(sequence)
        if len(sequence) < window_size or window_size == 0:
            result.append(f">{title}\n{sequence}")
            return result
        for i in range(len(sequence) - window_size + 1):
            start = i + 1
            end = i + window_size
            header = f">{title}_{start}_{end}"
            fragment = sequence[i:end]
            result.append(f"{header}\n{fragment}")
        return result

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if title is not None:
                fragments = process_sequence(title, sequence)
                for frag in fragments:
                    current_group.append(frag)
                    if len(current_group) == group_size:
                        fasta_groups.append("\n".join(current_group))
                        current_group = []
            title = line[1:].strip()
            sequence = ""
        else:
            sequence += line

    # 處理最後一條序列
    if title:
        fragments = process_sequence(title, sequence)
        for frag in fragments:
            current_group.append(frag)
            if len(current_group) == group_size:
                fasta_groups.append("\n".join(current_group))
                current_group = []

    if current_group:
        fasta_groups.append("\n".join(current_group))

    return fasta_groups, original_lengths
#run blastp
def run_blastp(fasta,db,evalue,matrix,gap):
    if gap == '0' :

        result = subprocess.run(
            [
                "blastp", 'blastp-short',
                "-query", "-",
                "-db", db,
                "-outfmt", "5",
                "-evalue", str(evalue),
                "-matrix", matrix,
                "-ungapped",
                "-comp_based_stats", "0",
                "-word_size", "2",
                "-out", "-"  # 輸出到 stdout
            ],
            input=fasta.encode(),  # 直接以 FASTA 字符串作為輸入
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
    else:


        result = subprocess.run(
            [
                "blastp",
                "-db", db,
                "-outfmt", "5",
                "-evalue", str(evalue),
                "-matrix", matrix,

                "-word_size", "2",
                "-out", "-"  # 輸出到 stdout
            ],
            input=fasta.encode(),  # 直接以 FASTA 字符串作為輸入
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
    blast_text = result.stdout.decode('utf-8')
    df = parse_blast_xml_to_dataframe(blast_text)
    return df
#xlm to csv
def parse_blast_xml_to_dataframe(xml_text):
    """將 BLAST XML 結果轉為 pandas DataFrame"""
    root = ET.fromstring(xml_text)
    records = []

    for iteration in root.find('BlastOutput_iterations').findall('Iteration'):
        query_id = iteration.findtext("Iteration_query-ID")
        query_def = iteration.findtext("Iteration_query-def")
        query_len = iteration.findtext("Iteration_query-len")

        hits = iteration.find("Iteration_hits")
        if hits is not None and hits.find("Hit") is not None:
            for hit in hits.findall("Hit"):
                hit_num = hit.findtext("Hit_num")
                hit_def = hit.findtext("Hit_def")
                hit_len = hit.findtext("Hit_len")

                hsps = hit.find("Hit_hsps")
                if hsps is not None:
                    for hsp in hsps.findall("Hsp"):
                        bit_score = hsp.findtext("Hsp_bit-score")
                        score = hsp.findtext("Hsp_score")
                        evalue = hsp.findtext("Hsp_evalue")
                        identity = hsp.findtext("Hsp_identity")
                        positive = hsp.findtext("Hsp_positive")
                        gaps = hsp.findtext("Hsp_gaps")
                        align_len = hsp.findtext("Hsp_align-len")

                        qseq = hsp.findtext("Hsp_qseq")
                        hseq = hsp.findtext("Hsp_hseq")
                        midline = hsp.findtext("Hsp_midline")
                        query_from = hsp.findtext("Hsp_query-from")
                        query_to = hsp.findtext("Hsp_query-to")
                        hit_from = hsp.findtext("Hsp_hit-from")
                        hit_to = hsp.findtext("Hsp_hit-to")
                        identity_percent = round(int(identity) / int(align_len) * 100, 1)
                        positive_percent = round(int(positive) / int(align_len) * 100, 1)
                        gaps_percent = round(int(gaps) / int(align_len) * 100, 1)
                        Query_MME_Len = int(query_to) - int(query_from) + 1
                        Hit_MME_Len=int(hit_to) - int(hit_from) + 1
                        records.append([
                            query_id, query_def, query_len,
                            hit_num, hit_def, hit_len,
                            bit_score, score, evalue,
                            identity, positive, gaps,identity_percent, positive_percent, gaps_percent, Query_MME_Len,
                            query_from, query_to, align_len,Hit_MME_Len, hit_from, hit_to,
                            qseq, hseq, midline
                        ])

    columns = [
        "Query_ID","Query_Def", "Query_Len",
        "Hit_Num", "Hit_Def", "Hit_Len",
        "Bit_Score", "Score", "Evalue",
        "Identity", "Positive", "Gaps","Identity%", "Positive%", "Gaps%", "Query_MME_Len",
        "Query_From", "Query_To", "Align_Len","Hit_MME_Len","Hit_From", "Hit_To",
        "Query_Seq", "Hit_Seq", "Midline"
    ]
    return pd.DataFrame(records, columns=columns)
##多工blastp
def process_one_group(fasta_chunk, blast_db, evalue, matrix, gap):
        try:
            blast_results = run_blastp(fasta_chunk, blast_db, evalue, matrix, gap)

            return blast_results
        except Exception as e:
            return f"[ERROR] {e}"
def run_all_groups_multicore(fasta_groups, blast_db, evalue, matrix, gap):
        all_results = []

        with ProcessPoolExecutor() as executor:
            futures = [
                executor.submit(process_one_group, chunk, blast_db, evalue, matrix, gap)
                for chunk in fasta_groups
            ]

            for i, future in enumerate(as_completed(futures), 1):
                try:
                    result = future.result()
                    print(f"✅ Group {i} finished.")
                    # print(result)
                    all_results.append(result)
                except Exception as e:
                    print(f"❌ Error in group {i}: {e}")
                    all_results.append(None)

            # 過濾掉 None 或錯誤訊息，只保留 DataFrame
            valid_results = [df for df in all_results if isinstance(df, pd.DataFrame)]

            # 合併成一張總表
            final_df = pd.concat(valid_results, ignore_index=True)
            # ✅ 如果結果為空就直接回傳
            if final_df.empty:
                print("⚠️ No results found, returning empty DataFrame.")
                return final_df
        return final_df
####query座標處理####
def adjust_query_coordinates(df, kmer):
    if kmer == 0:
        # 不做修正
        return df

    def extract_base(def_str):
        match = re.search(r"_(\d+)_(\d+)$", def_str)
        return int(match.group(1)) if match else 0

    def compute_start(row):
        base = extract_base(row["Query_Def"])
        return base + int(row["Query_From"]) - 1

    def compute_end(row):
        base = extract_base(row["Query_Def"])
        return base + int(row["Query_To"]) - 1

    df["Query_From"] = df.apply(compute_start, axis=1)
    df["Query_To"] = df.apply(compute_end, axis=1)

    return df
#根據 Query_Def 計算原始長度並取代 Query_Len 欄位
def extract_original_def(def_str):
    match = re.match(r"(.+?)_(\d+)_(\d+)$", def_str)
    return match.group(1) if match else def_str
#把 Query_Def 上的座標去掉
def strip_suffix_if_kmer(query_def, kmer):
        if int(kmer) > 0:
            # 匹配 _數字_數字 結尾
            return re.sub(r'_(\d+)_(\d+)$', '', query_def)
        return query_def
#去重複並把Evalue最小的留下
def dedup_stay_small_evalue(data):
    data["Evalue"] = pd.to_numeric(data["Evalue"], errors="coerce")  # 確保是數值
    final_df = (
        data.sort_values("Evalue")  # 先排序，Evalue小的會排前面
            .groupby(["Query_Def","Query_From" ,"Query_To","Hit_Def", "Hit_From", "Hit_To", "Hit_Seq"], as_index=False)
            .first()  # 每組只保留 Evalue 最小的那筆
        )
    return final_df
#二次篩選(長度範圍、identity、positive、gap)
def second_filter(blastp_data,min_length,max_length,identity,positive,gap_num):
            # 將 identity 從百分比轉成小數
            identity_threshold = identity * 0.01
            positive_threshold = positive * 0.01
            gap_threshold = gap_num * 0.01
            # 先把 Align_Len 轉成整數
            blastp_data['Align_Len'] = blastp_data['Align_Len'].astype(int)

            # 如果 Identity 也是字串，也一併轉成數字
            blastp_data['Identity'] = blastp_data['Identity'].astype(float)
            blastp_data['Positive'] = blastp_data['Positive'].astype(float)
            blastp_data['Gaps'] = blastp_data['Gaps'].astype(float)
            # 建立 mask，預設全 True
            mask = pd.Series(True, index=blastp_data.index)
            # 篩選最小長度
            if min_length > 0:
                mask &= blastp_data['Align_Len'] >= min_length

            # 篩選最大長度
            if max_length > 0:
                mask &= blastp_data['Align_Len'] <= max_length

            # 篩選 Identity/Align_Len 大於閾值
            if identity > 0:
                mask &= ((blastp_data['Identity'] / blastp_data['Align_Len']).round(2)) >= identity_threshold
            # 篩選 positive/Align_Len 大於閾值
            if positive > 0:
                mask &= ((blastp_data['Positive'] / blastp_data['Align_Len']).round(2)) >= positive_threshold
            # 篩選 Identity/Align_Len 大於閾值
            if gap_num > 0:
                mask &= ((blastp_data['Gaps'] / blastp_data['Align_Len']).round(2)) <= gap_threshold

            # 套用篩選
            blastp_data = blastp_data[mask]

            return blastp_data


#############################################################
#                     主 Pipeline
#############################################################
def blastp_pipeline_main(job_id):
    global CURRENT_JOB_DIR

    # === 專案 jobs 根目錄 ===
    project_root = Path(__file__).resolve().parent.parent
    out_dir = project_root / "jobs" / job_id
    job_dir = out_dir
    CURRENT_JOB_DIR = out_dir

    # === Step 0: 設定 log ===
    log_path = job_dir / "pipeline.log"

    class Tee:
        def __init__(self, *files):
            self.files = files
        def write(self, data):
            for f in self.files:
                f.write(data)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()

    log_file = open(log_path, "w", encoding="utf-8")
    sys.stdout = Tee(sys.stdout, log_file)
    sys.stderr = Tee(sys.stderr, log_file)

    success = False

    try:
        print("📝 Log 檔案建立成功：", log_path)
        print("🚀 Pipeline 開始執行...\n")

        # === 讀取設定檔 ===
        config_path = job_dir / "config.yaml"
        if not config_path.exists():
            raise FileNotFoundError(f"找不到設定檔: {config_path}")

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        logic = config.get("Matching_Logic")
        if logic == "Blastp":
            blast_db = Path(config["blastp_paths"]["human_db_path"]).resolve()  
            query_fasta_path = Path(config["blastp_paths"]["query_fasta_path"]).resolve()
            sequence_text = query_fasta_path.read_text(encoding="utf-8")
            human_location_path = Path(config["blastp_paths"]["human_protein_detail_path"]).resolve()
            IEDB_path = Path(config["blastp_paths"]["IEDB_file_path"]).resolve()

            evalue = config["blastp_parameters"]["evalue"]
            matrix = config["blastp_parameters"]["matrix"]
            gap = config["blastp_parameters"]["gap"]
            kmer = config["blastp_parameters"]["kmer_blastp"]
            min_length = config["blastp_parameters"]["min_length"]
            max_length = config["blastp_parameters"]["max_length"]
            identity = config["blastp_parameters"]["identity"]
            positive = config["blastp_parameters"]["positive"]
            gap_num = config["blastp_parameters"]["gap_num"]

        print("Step 1：切割 FASTA 序列中 ...")
        fasta_groups, original_lengths = sliding_windows_fasta(sequence_text, kmer, group_size=100)
        print(f"Step 1 完成：產生 {len(fasta_groups)} 組 FASTA 切片")

        print("Step 2：執行多核心 BLASTP ...")
        blastp_df = run_all_groups_multicore(fasta_groups, blast_db, evalue, matrix, gap)
        print("Step 2 完成")

        print("Step 3：調整座標與 Query_Len ...")
        blastp_df_coordinate = adjust_query_coordinates(blastp_df, kmer)
        blastp_df_coordinate["Query_Len"] = blastp_df_coordinate["Query_Def"].apply(
            lambda qdef: original_lengths.get(extract_original_def(qdef), None)
        )
        blastp_df_coordinate["Query_Def"] = blastp_df_coordinate["Query_Def"].apply(lambda x: strip_suffix_if_kmer(x, kmer))
        print("Step 3 完成")

        print("Step 4：去重複並保留最小 Evalue ...")
        blastp_df_coordinate_dedup = dedup_stay_small_evalue(blastp_df_coordinate)
        print("Step 4 完成")

        # 陳昱丞修改部分
        filtered_df = second_filter(
            blastp_df_coordinate_dedup,
            min_length, max_length, identity, positive, gap_num
        )
        
        filtered_df = filtered_df.rename(columns={
            "Query_Def": "query_protein_name",
            "Query_From": "MME(query)_start",
            "Query_To": "MME(query)_end",
            "Hit_Def": "hit_human_protein_name",
            "Hit_From": "MME(hit)_start",
            "Hit_To": "MME(hit)_end",
            "Hit_Seq": "MME(hit)",
            "Query_Len": "query_protein_length",
            "Hit_Len": "hit_human_protein_length",
            "Query_MME_Len": "length_of_MME(query)",
            "Hit_MME_Len": "length_of_MME(hit)",
            "Query_Seq": "MME(query)",
        })
        filtered_df.drop(columns=["Query_ID","Hit_Num"], inplace=True)

        core_columns = [
            "MME(query)",
            "MME(hit)",
            "query_protein_name",
            "query_protein_length",
            "length_of_MME(query)",
            "MME(query)_start",
            "MME(query)_end",
            "hit_human_protein_name",
            "hit_human_protein_length",
            "length_of_MME(hit)",
            "MME(hit)_start",
            "MME(hit)_end",
        ]
        for col in core_columns:
            if col not in filtered_df.columns:
                filtered_df[col] = None   # 或 np.nan
        other_columns = [
            col for col in filtered_df.columns
            if col not in core_columns
        ]
        filtered_df = filtered_df[core_columns + other_columns]


        print("perfect match分析完成")
        IEDB_match_result = IEDB_analyze_matches(filtered_df, IEDB_path, out_dir)
        calc_reference_table(IEDB_match_result, human_location_path, out_dir)
        calc_epitope_table(IEDB_match_result, out_dir)
        calc_query_table(IEDB_match_result, out_dir)
        print("Pipeline 完成！")
        return filtered_df

    except Exception as e:
        error_message = traceback.format_exc()
        print("❌ Pipeline 發生錯誤：\n", error_message)

    finally:
        # 恢復 stdout/stderr
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        log_file.close()
        