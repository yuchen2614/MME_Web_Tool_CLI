import pandas as pd
import time
from collections import defaultdict
from io import StringIO
import yaml
import sys
import traceback
from pathlib import Path
from MME_Web_Tool_CLI.analysis_pipeline_code.shared_analysis_steps import \
    IEDB_analyze_matches, calc_reference_table, calc_epitope_table, calc_query_table

# perfect_match分析
def find_perfect_matches_df(human_fasta: Path, virus_fasta: Path,k: int, job_dir: Path) -> pd.DataFrame:

    def read_fasta_headers_and_seqs(fasta_input):
        """
        支援：
        - Path 物件
        - 檔案路徑字串
        - 直接傳入 FASTA 字串內容
        """
        headers, seqs = [], []
        buf, head = [], None

        # 轉成 Path
        if isinstance(fasta_input, (str, Path)) and Path(fasta_input).exists():
            f = Path(fasta_input).open("r", encoding="utf-8", errors="ignore")
        elif isinstance(fasta_input, str):
            f = StringIO(fasta_input)
        else:
            raise TypeError("fasta_input 必須是 Path 或檔案路徑字串")

        with f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                if s.startswith(">"):
                    if head is not None:
                        seqs.append("".join(buf))
                        buf = []
                    head = s[1:]
                    headers.append(head)
                else:
                    buf.append(s)

            if head is not None:
                seqs.append("".join(buf))

        return headers, seqs

    # 開始分析
    start_time = time.perf_counter()

    human_headers, human_seqs = read_fasta_headers_and_seqs(human_fasta)
    virus_headers, virus_seqs = read_fasta_headers_and_seqs(virus_fasta)

    human_len_by_header = {h: len(s) for h, s in zip(human_headers, human_seqs)}
    human_seq_by_header = {h: s for h, s in zip(human_headers, human_seqs)}

    # 先收集病毒所有 k-mer
    virus_kset = set()
    for vs in virus_seqs:
        for b in range(0, max(0, len(vs) - k + 1)):
            virus_kset.add(vs[b:b+k])

    # 建人類 k-mer 索引
    human_index = defaultdict(list)
    for hh, hs in zip(human_headers, human_seqs):
        for a in range(0, max(0, len(hs) - k + 1)):
            kmer = hs[a:a+k]
            if kmer in virus_kset:
                human_index[kmer].append((hh, a))

    rows = []

    for vhead, vseq in zip(virus_headers, virus_seqs):
        Lq = len(vseq)
        active = defaultdict(list)

        for b in range(0, max(0, Lq - k + 1)):
            kmer = vseq[b:b+k]
            touched_keys = set()

            for hhead, a in human_index.get(kmer, ()):
                key = (hhead, a - b)
                lst = active.get(key, [])
                extended = False

                for seg in lst:
                    if seg["next_b"] == b and seg["next_a"] == a:
                        seg["len"] += 1
                        seg["next_b"] += 1
                        seg["next_a"] += 1
                        extended = True
                        break

                if not extended:
                    lst.append({"qb": b, "hb": a, "next_b": b+1, "next_a": a+1, "len": k})
                    active[key] = lst

                touched_keys.add(key)

            # flush 沒有延續的段落
            for key in list(active.keys()):
                if key not in touched_keys:
                    hhead = key[0]
                    hseq  = human_seq_by_header[hhead]
                    hlen  = human_len_by_header[hhead]
                    for seg in active.pop(key):
                        q0, h0, ln = seg["qb"], seg["hb"], seg["len"]
                        rows.append({
                            'MME(query)': vseq[q0:q0+ln],
                            'MME(hit)':   hseq[h0:h0+ln],
                            'query_protein_name': vhead,
                            'query_protein_length': Lq,
                            'length_of_MME(query)': ln,
                            'MME(query)_start': q0 + 1,
                            'MME(query)_end':   q0 + ln,
                            'hit_human_protein_name': hhead,
                            'hit_human_protein_length': hlen,
                            'length_of_MME(hit)': ln,
                            'MME(hit)_start': h0 + 1,
                            'MME(hit)_end':   h0 + ln,
                        })

        # 掃描結束，把剩餘段落輸出
        for key, lst in active.items():
            hhead = key[0]
            hseq  = human_seq_by_header[hhead]
            hlen  = human_len_by_header[hhead]
            for seg in lst:
                q0, h0, ln = seg["qb"], seg["hb"], seg["len"]
                rows.append({
                    'MME(query)': vseq[q0:q0+ln],
                    'MME(hit)':   hseq[h0:h0+ln],
                    'query_protein_name': vhead,
                    'query_protein_length': Lq,
                    'length_of_MME(query)': ln,
                    'MME(query)_start': q0 + 1,
                    'MME(query)_end':   q0 + ln,
                    'hit_human_protein_name': hhead,
                    'hit_human_protein_length': hlen,
                    'length_of_MME(hit)': ln,
                    'MME(hit)_start': h0 + 1,
                    'MME(hit)_end':   h0 + ln,
                })
    # ===== 新增 Epitope_unique_in_reference =====
    df = pd.DataFrame(rows)
    print(f"perfect match分析完成，用時 {time.perf_counter() - start_time:.2f} 秒。\n")
    return df


# 完整分析主程式
def perfect_match_pipeline_main(job_id):

    project_root = Path(__file__).resolve().parent.parent
    job_dir = project_root / "jobs" / job_id

    # === 設定 pipeline.log ===
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

    try:
        print("📝 Log 檔案建立成功：", log_path)
        print("🚀 Pipeline 開始執行...\n")

        # === 讀取 job_id 資料夾的設定檔 config.yaml ===
        config_path = job_dir / "config.yaml"
        if not config_path.exists():
            raise FileNotFoundError(f"找不到設定檔: {config_path}")

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        logic = config.get("Matching_Logic")
        if logic == "Perfect_Match":

            human_fasta_path = Path(config["perfect_match_paths"]["human_fasta_path"]).resolve()
            query_fasta_path = Path(config["perfect_match_paths"]["query_fasta_path"]).resolve()
            human_protein_detail_path = Path(config["perfect_match_paths"]["human_protein_detail_path"]).resolve()
            IEDB_file_path = Path(config["perfect_match_paths"]["IEDB_file_path"]).resolve()

            kmer_perfect_match = config["perfect_match_parameters"]["kmer_perfect_match"]

        # === Pipeline Steps ===
        print("🚀 Step 1: Perfect Match 分析中...")
        perfect_match_result_df = find_perfect_matches_df(human_fasta_path, query_fasta_path, kmer_perfect_match, job_dir)
        if perfect_match_result_df.empty:
            raise ValueError("⚠️ 找不到任何 perfect match 結果，Pipeline 結束。")

        print("🚀 Step 2: IEDB 分析中...")
        perfect_match_result_df = IEDB_analyze_matches(perfect_match_result_df, IEDB_file_path, job_dir)

        print("🚀 Step 3: Reference 統計中...")
        calc_reference_table(perfect_match_result_df, human_protein_detail_path, job_dir)

        print("🚀 Step 4: Epitope 統計中...")
        calc_epitope_table(perfect_match_result_df, job_dir)

        print("🚀 Step 5: Query 統計中...")
        calc_query_table(perfect_match_result_df, job_dir)

        print(f"✅ Pipeline 完成！結果輸出：{job_dir}")
        print(f"📁 Log 檔案位置：{log_path}")

    except Exception as e:
        error_message = traceback.format_exc()
        print("❌ Pipeline 發生錯誤：\n", error_message)

    finally:
        # 恢復 stdout/stderr
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        log_file.close()













    