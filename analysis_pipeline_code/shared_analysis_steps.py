import pandas as pd
import time
import numpy as np
import ahocorasick
from pathlib import Path

# IEDB分析
def IEDB_analyze_matches(match_df: pd.DataFrame, IEDB_file_path: Path, job_dir: Path):
    FAST_MODE_THRESHOLD = 10_000
    USE_FAST_MODE = len(match_df) < FAST_MODE_THRESHOLD
    print(f"IEDB: rows={len(match_df)}, 使用 {'高速 merge 模式' if USE_FAST_MODE else '低記憶體分組模式'}")

    start_time = time.perf_counter()

    # --- 確保為 Path 物件並轉為絕對路徑 ---
    IEDB_file_path = Path(IEDB_file_path).resolve()
    # --- 讀取 IEDB CSV ---
    IEDB_human_df = pd.read_csv(IEDB_file_path)
    df = match_df.copy()
    df["hit_human_protein_id"] = df["hit_human_protein_name"].str.split("|", expand=True)[1]


    # ================================================================
    # PART 1 — substring count（Aho-Corasick）
    # ================================================================
    iedb_seqs = sorted(IEDB_human_df["Name"].dropna().unique(), key=len)

    # 🔸 把 hit 跟 query 的 MME 都合併，避免建兩棵樹
    unique_mmes_hit = df["MME(hit)"].dropna().unique()
    unique_mmes_query = df["MME(query)"].dropna().unique()

    # 建立全集避免漏掉某邊的 MME
    all_unique_mmes = set(unique_mmes_hit) | set(unique_mmes_query)

    # --- Build Aho-Corasick automaton ---
    A = ahocorasick.Automaton()
    for seq in all_unique_mmes:
        A.add_word(seq, seq)
    A.make_automaton()

    # --- Count matches ---
    counts = dict.fromkeys(all_unique_mmes, 0)

    for epitope in iedb_seqs:
        for _, match in A.iter(epitope):
            counts[match] += 1

    # --- 分別寫入兩欄 ---
    df["IEDB_query_epitope_substring_count"] = df["MME(query)"].map(counts).fillna(0).astype(int)
    df["IEDB_human_epitope_substring_count"] = df["MME(hit)"].map(counts).fillna(0).astype(int)

    # ================================================================
    # 準備 IEDB 資料
    # ================================================================
    IEDB_human = IEDB_human_df.rename(columns={
        "UniProt_ID": "hit_human_protein_id",
        "Starting Position": "IEDB_human_start",
        "Ending Position": "IEDB_human_end"
    }).copy()

    IEDB_human["IEDB_human_start"] = IEDB_human["IEDB_human_start"].fillna(0).astype(int)
    IEDB_human["IEDB_human_end"] = IEDB_human["IEDB_human_end"].fillna(0).astype(int)
    IEDB_human["IRI_base"] = IEDB_human["IEDB IRI"].astype(str).apply(lambda x: x.split("/")[-1])

    # ================================================================
    # MODE A — 高速 merge（小 df）
    # ================================================================
    if USE_FAST_MODE:
        df_reset = df.reset_index().rename(columns={"index": "df_idx"})
        merged = df_reset.merge(
            IEDB_human[["hit_human_protein_id", "IEDB_human_start", "IEDB_human_end", "IRI_base"]],
            on="hit_human_protein_id",
            how="left"
        )

        MS = merged["MME(hit)_start"].fillna(0).astype(int)
        ME = merged["MME(hit)_end"].fillna(0).astype(int)
        IS = merged["IEDB_human_start"].fillna(0).astype(int)
        IE = merged["IEDB_human_end"].fillna(0).astype(int)

        mask_fully = (MS >= IS) & (ME <= IE)
        mask_partial = (ME >= IS) & (MS <= IE) & (~mask_fully)

        merged["fully"] = mask_fully
        merged["partial"] = mask_partial

        # fully 聚合
        agg_fully = merged.groupby("df_idx").agg({
            "fully": "sum",
            "IRI_base": lambda x: ", ".join(
                x[merged.loc[x.index, "fully"]]
            ) if any(merged.loc[x.index, "fully"]) else "No_evidence"
        }).rename(columns={
            "fully": "IEDB_human_positional_fully_contained",
            "IRI_base": "IEDB_fully_contained_IRI"
        })

        # partial 聚合
        agg_partial = merged.groupby("df_idx").agg({
            "partial": "sum",
            "IRI_base": lambda x: ", ".join(
                x[merged.loc[x.index, "partial"]]
            ) if any(merged.loc[x.index, "partial"]) else "No_evidence"
        }).rename(columns={
            "partial": "IEDB_human_positional_partial_overlap",
            "IRI_base": "IEDB_partial_overlap_IRI"
        })

        df = df_reset.merge(agg_fully, left_on="df_idx", right_index=True, how="left")
        df = df.merge(agg_partial, left_on="df_idx", right_index=True, how="left")

        df.drop(columns=["df_idx"], inplace=True)

    # ================================================================
    # MODE B — 低記憶體（大 df）
    # ================================================================
    else:
        # ================================================================
        # 初始化欄位
        # ================================================================
        df["IEDB_human_positional_fully_contained"] = 0
        df["IEDB_human_positional_partial_overlap"] = 0
        df["IEDB_fully_contained_IRI"] = "No_evidence"
        df["IEDB_partial_overlap_IRI"] = "No_evidence"

        for protein_id, sub_df in df.groupby("hit_human_protein_id", group_keys=False):

            sub_iedb = IEDB_human[IEDB_human["hit_human_protein_id"] == protein_id]
            if sub_iedb.empty:
                continue

            I_st = sub_iedb["IEDB_human_start"].astype(int).to_numpy()
            I_en = sub_iedb["IEDB_human_end"].astype(int).to_numpy()
            I_iri = sub_iedb["IRI_base"].astype(int).to_numpy()

            M_st = sub_df["MME(hit)_start"].astype(int).to_numpy()
            M_en = sub_df["MME(hit)_end"].astype(int).to_numpy()
            idxs = sub_df.index

            fully_cnt, partial_cnt = [], []
            fully_iri_list, partial_iri_list = [], []

            for ms, me in zip(M_st, M_en):
                fully_mask = (ms >= I_st) & (me <= I_en)
                partial_mask = (me >= I_st) & (ms <= I_en) & (~fully_mask)

                if fully_mask.any():
                    fully_cnt.append(int(fully_mask.sum()))
                    fully_iri_list.append(", ".join(I_iri[fully_mask]))
                else:
                    fully_cnt.append(0)
                    fully_iri_list.append("No_evidence")

                if partial_mask.any():
                    partial_cnt.append(int(partial_mask.sum()))
                    partial_iri_list.append(", ".join(I_iri[partial_mask]))
                else:
                    partial_cnt.append(0)
                    partial_iri_list.append("No_evidence")

            df.loc[idxs, "IEDB_human_positional_fully_contained"] = fully_cnt
            df.loc[idxs, "IEDB_fully_contained_IRI"] = fully_iri_list
            df.loc[idxs, "IEDB_human_positional_partial_overlap"] = partial_cnt
            df.loc[idxs, "IEDB_partial_overlap_IRI"] = partial_iri_list

    # ================================================================
    # human protein epitope count
    # ================================================================
    human_count_dict = IEDB_human["hit_human_protein_id"].value_counts().to_dict()
    df["IEDB_human_protein_data_count"] = df["hit_human_protein_id"].map(human_count_dict).fillna(0).astype(int)

    # evidence 欄位
    df["IEDB_evidence"] = np.where(
        df["IEDB_human_positional_fully_contained"] > 0, "Strong_evidence",
        np.where(df["IEDB_human_positional_partial_overlap"] > 0, "Weak_evidence", "No_evidence")
    )

    # ================================================================
    # 欄位順序
    # ================================================================
    tail_cols = [
        "IEDB_human_protein_data_count",
        "IEDB_human_positional_fully_contained",
        "IEDB_human_positional_partial_overlap",
        "IEDB_fully_contained_IRI",
        "IEDB_partial_overlap_IRI",
        "IEDB_evidence"
    ]
    front_cols = [c for c in df.columns if c not in tail_cols]
    df = df[front_cols + tail_cols]

    df.to_csv(f"{job_dir}/perfect_match.csv", index=False)
    print(f"IEDB分析完成，用時 {time.perf_counter() - start_time:.2f} 秒。\n")
    return df

# reference統計表
def calc_reference_table(perfect_match_result_df: pd.DataFrame, human_protein_detail_path: Path, job_dir: Path):
    start_time = time.perf_counter()

    # --- 確保為 Path 物件並轉為絕對路徑 ---
    human_protein_detail_path = Path(human_protein_detail_path).resolve()
    job_dir = Path(job_dir)

    # --- 讀取 Human Protein Detail CSV ---
    human_protein_detail_df = pd.read_csv(human_protein_detail_path)

    df = perfect_match_result_df.copy()

    # 計算統計資訊
    agg_df = perfect_match_result_df.groupby('hit_human_protein_id').agg(
        epitope_count=('MME(query)', 'size'),
        query_protein_name_kind=('query_protein_name', 'nunique'),
        Strong_IEDB_evidence_MME=('IEDB_evidence', lambda x: (x=='Strong_evidence').sum()),
        Weak_IEDB_evidence_MME=('IEDB_evidence', lambda x: (x=='Weak_evidence').sum())
    ).reset_index()

    # 將 human_protein_detail_df 只保留必要欄位
    human_protein_detail_df = human_protein_detail_df[['UniProt_protein', 'Gene_HGNC', 'Gene_description','Main_location']]

    # merge 結合
    result_reference = agg_df.merge(
        human_protein_detail_df,
        left_on='hit_human_protein_id',
        right_on='UniProt_protein',
        how='left'
    )

    # 空值補 'undefined'
    result_reference['Gene_HGNC'] = result_reference['Gene_HGNC'].fillna('No data')
    result_reference['Gene_description'] = result_reference['Gene_description'].fillna('No data')
    result_reference['Main_location'] = result_reference['Main_location'].fillna('No data')
    result_reference = result_reference.rename(columns={"Gene_description": "Protein_name",})
    # 重新調整欄位順序：hit_human_protein_id, Gene_HGNC, Gene_description, epitope_count, ...
    cols = ['hit_human_protein_id', 'Protein_name', 'Gene_HGNC', 'Main_location',
            'epitope_count', 'query_protein_name_kind', 
            'Strong_IEDB_evidence_MME', 'Weak_IEDB_evidence_MME']
    result_reference = result_reference[cols]

    # 輸出 CSV
    result_reference.to_csv(f'{job_dir}/reference.csv', index=False)

    print(f"Reference分析完成，用時 {time.perf_counter() - start_time:.2f} 秒。\n")

# epitope統計表
def calc_epitope_table(perfect_match_result_df, job_dir):
    start_time = time.perf_counter()

    # groupby MME(hit)
    agg_df = perfect_match_result_df.groupby("MME(hit)").agg(
        epitope_length=("MME(hit)", lambda x: len(x.iloc[0])),
        epitope_count=("MME(hit)", "size"),
        hit_human_protein_id_kind=("hit_human_protein_id", "nunique"),
        query_protein_name_kind=("query_protein_name", "nunique"),
        IEDB_human_epitope_substring_count=("IEDB_human_epitope_substring_count", "first")
    ).reset_index()
    agg_df = agg_df.rename(columns={"MME(hit)": "Epitope"})
    agg_df.to_csv(f"{job_dir}/epitope.csv", index=False)
    # groupby MME(query)
    agg_df = perfect_match_result_df.groupby("MME(query)").agg(
        epitope_length=("MME(query)", lambda x: len(x.iloc[0])),
        epitope_count=("MME(query)", "size"),
        hit_human_protein_id_kind=("hit_human_protein_id", "nunique"),
        query_protein_name_kind=("query_protein_name", "nunique"),
    ).reset_index()
    agg_df = agg_df.rename(columns={"MME(query)": "Epitope"})
    agg_df.to_csv(f"{job_dir}/epitope_query.csv", index=False)
    print(f"Epitope分析完成，用時 {time.perf_counter() - start_time:.2f} 秒。\n")

# query統計表
def calc_query_table(perfect_match_result_df, job_dir):
    start_time = time.perf_counter()

    # groupby query_protein_name，一次計算所有統計值
    agg_df = perfect_match_result_df.groupby("query_protein_name").agg(
        hit_human_protein_id_kind=("hit_human_protein_id", "nunique"),
        epitope_count=("MME(hit)", "size")
    ).reset_index()

    # 輸出 CSV
    agg_df.to_csv(f"{job_dir}/query.csv", index=False)

    print(f"Query分析完成，用時 {time.perf_counter() - start_time:.2f} 秒。\n")
