#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_enso_cross_index_skill_persistence_html.py

Compara DUAS curvas climáticas, permitindo qualquer combinação entre ONI e RONI:

    referência ONI  vs comparação ONI
    referência ONI  vs comparação RONI
    referência RONI vs comparação ONI
    referência RONI vs comparação RONI

O script mantém a robustez de leitura da versão anterior e adiciona:

    1. Opção para salvar ou não salvar CSVs diagnósticos.
    2. Métricas categóricas de skill EN/LN/N.
    3. Índices de persistência ENSO para referência e comparação.
    4. Diferença comparativa dos índices de persistência entre as duas séries.
    5. HTML interativo com:
        - série temporal + pontos de divergência no topo;
        - matriz de confusão EN/N/LN;
        - tabela de métricas de skill;
        - análise de persistência e duração dos regimes ENSO.

Classificação ENSO implementada:
    EN = El Niño  : sequência contínua de >= 5 valores consecutivos > +0.5 °C
    LN = La Niña  : sequência contínua de >= 5 valores consecutivos < -0.5 °C
    N  = Neutro   : todos os demais casos

Ponto metodológico importante:
    Quando uma sequência tem 5 ou mais valores consecutivos acima/abaixo do limiar,
    TODOS os trimestres dessa sequência são classificados como EN/LN, não apenas
    o quinto trimestre.

Requisitos:
    pip install pandas numpy plotly
"""

# ============================================================
# 1. CONFIGURAÇÕES DO USUÁRIO
# ============================================================

# Diretório da série de referência.
REFERENCE_DIRECTORY = r"./resultados/ersstv5"

# Diretório da única série a ser comparada à referência.
COMPARISON_DIRECTORY = r"./resultados/ersstv5"

# Índice usado na referência: "ONI" ou "RONI".
REFERENCE_INDEX = "ONI"

# Índice usado na comparação: "ONI" ou "RONI".
COMPARISON_INDEX = "RONI"

# Rótulos opcionais. Se None, usa "ÍNDICE (nome_da_pasta)".
REFERENCE_LABEL = None
COMPARISON_LABEL = None

# Critérios ENSO.
WARM_THRESHOLD = 0.5
COLD_THRESHOLD = -0.5
MIN_CONSECUTIVE_SEASONS = 5

# ------------------------------------------------------------
# Saídas
# ------------------------------------------------------------
OUTPUT_HTML = "ENSO_cross_index_skill_persistence.html"

# Opção solicitada: reativar CSVs quando desejado.
# False -> gera apenas o HTML.
# True  -> gera HTML + CSVs diagnósticos.
SAVE_CSV = False

OUTPUT_ALIGNED_DIAGNOSTIC_CSV = "ENSO_cross_index_aligned_diagnostic.csv"
OUTPUT_CONFUSION_MATRIX_CSV = "ENSO_cross_index_confusion_matrix.csv"
OUTPUT_SKILL_METRICS_CSV = "ENSO_cross_index_skill_metrics.csv"
OUTPUT_EVENT_BLOCKS_CSV = "ENSO_cross_index_event_blocks.csv"
OUTPUT_PERSISTENCE_SUMMARY_CSV = "ENSO_cross_index_persistence_summary.csv"
OUTPUT_PERSISTENCE_COMPARISON_CSV = "ENSO_cross_index_persistence_comparison.csv"

# ------------------------------------------------------------
# Visual
# ------------------------------------------------------------
HTML_TITLE = "Diagnóstico ENSO: comparação entre ONI/RONI, skill e persistência"
PLOT_TEMPLATE = "plotly_white"
FIGURE_WIDTH = 1550
TIME_SERIES_HEIGHT = 780
CONFUSION_HEIGHT = 560
SKILL_HEIGHT = 520
PERSISTENCE_HEIGHT = 760

REFERENCE_LINE_COLOR = "black"
COMPARISON_LINE_COLOR = "darkorange"
REFERENCE_LINE_WIDTH = 2.3
COMPARISON_LINE_WIDTH = 1.9

SHOW_THRESHOLDS = True
SHOW_ZERO_LINE = True
TOP_BAND_FRACTION = 0.32
DIVERGENCE_MARKER_SIZE = 12
DECIMALS = 3

# Ordem das classes no diagnóstico.
CLASS_ORDER = ["EN", "N", "LN"]

# Cores dos tipos de divergência.
DIVERGENCE_COLORS = {
    "REF_EN__NEW_N": "#ff7f0e",      # laranja
    # "REF_EN__NEW_LN": "#9467bd",     # roxo
    "REF_N__NEW_EN": "#d62728",      # vermelho
    "REF_N__NEW_LN": "#1f77b4",      # azul
    "REF_LN__NEW_N": "#17becf",      # ciano
    # "REF_LN__NEW_EN": "#e377c2",     # magenta
}

DIVERGENCE_LABELS = {
    "REF_EN__NEW_N": "Referência EN / comparação N",
    # "REF_EN__NEW_LN": "Referência EN / comparação LN",
    "REF_N__NEW_EN": "Referência N / comparação EN",
    "REF_N__NEW_LN": "Referência N / comparação LN",
    "REF_LN__NEW_N": "Referência LN / comparação N",
    # "REF_LN__NEW_EN": "Referência LN / comparação EN",
}

DIVERGENCE_ORDER = [
    "REF_EN__NEW_N",
    # "REF_EN__NEW_LN",
    "REF_N__NEW_EN",
    "REF_N__NEW_LN",
    "REF_LN__NEW_N",
    # "REF_LN__NEW_EN",
]

FILE_MAP = {
    "ONI": "ONI_timeseries.csv",
    "RONI": "RONI_timeseries.csv",
}

VALUE_COL_MAP = {
    "ONI": "ONI",
    "RONI": "RONI",
}

# ============================================================
# 2. IMPORTS
# ============================================================

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import to_html


# ============================================================
# 3. VALIDAÇÃO E LEITURA ROBUSTA
# ============================================================

def validate_config() -> None:
    """Valida configurações principais."""
    for name, value in [("REFERENCE_INDEX", REFERENCE_INDEX), ("COMPARISON_INDEX", COMPARISON_INDEX)]:
        if value not in FILE_MAP:
            raise ValueError(f"{name} deve ser 'ONI' ou 'RONI'. Valor recebido: {value}")

    ref_dir = Path(REFERENCE_DIRECTORY).expanduser().resolve()
    cmp_dir = Path(COMPARISON_DIRECTORY).expanduser().resolve()

    if not ref_dir.exists() or not ref_dir.is_dir():
        raise FileNotFoundError(f"Diretório de referência não encontrado: {ref_dir}")

    if not cmp_dir.exists() or not cmp_dir.is_dir():
        raise FileNotFoundError(f"Diretório de comparação não encontrado: {cmp_dir}")


def make_label(directory: str, index_name: str, explicit_label: Optional[str]) -> str:
    """Define rótulo da série."""
    if explicit_label:
        return explicit_label
    p = Path(directory).expanduser().resolve()
    folder_name = p.name or str(p)
    return f"{index_name} ({folder_name})"


def read_index_file(directory: str, index_name: str) -> pd.DataFrame:
    """
    Leitura robusta de ONI_timeseries.csv ou RONI_timeseries.csv.

    Aceita:
        - coluna time; ou
        - colunas year e month.

    Retorna sempre:
        time, year, month, season, value
    """
    directory = Path(directory).expanduser().resolve()
    file_path = directory / FILE_MAP[index_name]
    value_col = VALUE_COL_MAP[index_name]

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

    df = pd.read_csv(file_path)
    df.columns = [str(c).strip() for c in df.columns]

    if value_col not in df.columns:
        raise KeyError(
            f"Coluna '{value_col}' não encontrada em {file_path}. "
            f"Colunas disponíveis: {list(df.columns)}"
        )

    # Tempo.
    if "time" in df.columns:
        df["time"] = pd.to_datetime(df["time"], errors="coerce")
    elif {"year", "month"}.issubset(df.columns):
        df["time"] = pd.to_datetime(
            dict(
                year=pd.to_numeric(df["year"], errors="coerce"),
                month=pd.to_numeric(df["month"], errors="coerce"),
                day=15,
            ),
            errors="coerce",
        )
    else:
        raise KeyError(f"O arquivo {file_path} precisa conter 'time' ou as colunas 'year' e 'month'.")

    if df["time"].isna().all():
        raise ValueError(f"Não foi possível interpretar a coluna temporal em {file_path}.")

    # year/month.
    if "year" not in df.columns:
        df["year"] = df["time"].dt.year
    else:
        df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")

    if "month" not in df.columns:
        df["month"] = df["time"].dt.month
    else:
        df["month"] = pd.to_numeric(df["month"], errors="coerce").astype("Int64")

    # season.
    if "season" not in df.columns:
        season_by_month = {
            1: "DJF", 2: "JFM", 3: "FMA", 4: "MAM",
            5: "AMJ", 6: "MJJ", 7: "JJA", 8: "JAS",
            9: "ASO", 10: "SON", 11: "OND", 12: "NDJ",
        }
        df["season"] = df["month"].astype("Int64").map(season_by_month)
    else:
        df["season"] = df["season"].astype(str).str.strip()

    df["value"] = pd.to_numeric(df[value_col], errors="coerce")

    df = df[["time", "year", "month", "season", "value"]].copy()
    df = df.dropna(subset=["time", "year", "month"])
    df["year"] = df["year"].astype(int)
    df["month"] = df["month"].astype(int)
    df = df.sort_values(["year", "month"]).reset_index(drop=True)

    return df


# ============================================================
# 4. CLASSIFICAÇÃO ENSO POR SEQUÊNCIAS COMPLETAS
# ============================================================

def classify_enso_full_runs(values: np.ndarray) -> List[str]:
    """
    Classifica cada trimestre como EN, LN ou N.

    Sequências contínuas completas com >= MIN_CONSECUTIVE_SEASONS recebem
    EN ou LN em todos os seus elementos.
    """
    values = np.asarray(values, dtype=float)
    labels = np.array(["N"] * len(values), dtype=object)

    warm_mask = np.isfinite(values) & (values > WARM_THRESHOLD)
    labels = apply_run_classification(labels, warm_mask, "EN")

    cold_mask = np.isfinite(values) & (values < COLD_THRESHOLD)
    labels = apply_run_classification(labels, cold_mask, "LN")

    return labels.tolist()


def apply_run_classification(labels: np.ndarray, mask: np.ndarray, class_label: str) -> np.ndarray:
    """Aplica class_label a sequências True com tamanho >= mínimo."""
    out = labels.copy()
    n = len(mask)
    i = 0

    while i < n:
        if not mask[i]:
            i += 1
            continue

        start = i
        while i < n and mask[i]:
            i += 1
        end = i - 1

        if (end - start + 1) >= MIN_CONSECUTIVE_SEASONS:
            out[start : end + 1] = class_label

    return out


def add_classification(df: pd.DataFrame) -> pd.DataFrame:
    """Adiciona classificação ENSO."""
    out = df.copy()
    out["enso_class"] = classify_enso_full_runs(out["value"].to_numpy())
    return out


# ============================================================
# 5. ALINHAMENTO, EVENTOS E PERSISTÊNCIA
# ============================================================

def align_reference_and_comparison(ref_df: pd.DataFrame, cmp_df: pd.DataFrame) -> pd.DataFrame:
    """Alinha referência e comparação por year/month/season."""
    ref = ref_df.rename(columns={"value": "value_ref", "enso_class": "class_ref"})
    cmp_ = cmp_df.rename(columns={"value": "value_cmp", "enso_class": "class_cmp"})

    aligned = pd.merge(
        ref[["time", "year", "month", "season", "value_ref", "class_ref"]],
        cmp_[["year", "month", "season", "value_cmp", "class_cmp"]],
        on=["year", "month", "season"],
        how="inner",
    )

    aligned["is_divergent"] = aligned["class_ref"] != aligned["class_cmp"]
    aligned["divergence_type"] = aligned.apply(classify_divergence_type, axis=1)
    return aligned


def classify_divergence_type(row: pd.Series) -> Optional[str]:
    """Código da divergência ou None."""
    if row["class_ref"] == row["class_cmp"]:
        return None
    return f"REF_{row['class_ref']}__NEW_{row['class_cmp']}"


def extract_class_blocks(df: pd.DataFrame, label: str, index_name: str) -> pd.DataFrame:
    """Extrai blocos contínuos EN/N/LN para análise de persistência."""
    if df.empty:
        return pd.DataFrame()

    rows = []
    current = df.iloc[0]["enso_class"]
    start_idx = 0

    for i in range(1, len(df)):
        cls = df.iloc[i]["enso_class"]
        if cls != current:
            rows.append(make_block_record(df, start_idx, i - 1, current, label, index_name))
            current = cls
            start_idx = i

    rows.append(make_block_record(df, start_idx, len(df) - 1, current, label, index_name))
    return pd.DataFrame(rows)


def make_block_record(df: pd.DataFrame, start_idx: int, end_idx: int, cls: str, label: str, index_name: str) -> Dict[str, object]:
    """Registro de bloco contínuo."""
    start = df.iloc[start_idx]
    end = df.iloc[end_idx]
    block = df.iloc[start_idx : end_idx + 1]
    return {
        "series": label,
        "index": index_name,
        "class": cls,
        "start_time": start["time"],
        "end_time": end["time"],
        "start_year": int(start["year"]),
        "start_month": int(start["month"]),
        "start_season": start["season"],
        "end_year": int(end["year"]),
        "end_month": int(end["month"]),
        "end_season": end["season"],
        "duration_seasons": int(end_idx - start_idx + 1),
        "mean_value": float(block["value"].mean(skipna=True)),
        "min_value": float(block["value"].min(skipna=True)),
        "max_value": float(block["value"].max(skipna=True)),
    }


def persistence_summary(blocks: pd.DataFrame) -> pd.DataFrame:
    """Resumo de persistência/duração por série, índice e classe."""
    if blocks.empty:
        return pd.DataFrame()

    summary = (
        blocks
        .groupby(["series", "index", "class"], as_index=False)
        .agg(
            n_blocks=("duration_seasons", "count"),
            total_seasons=("duration_seasons", "sum"),
            mean_duration=("duration_seasons", "mean"),
            median_duration=("duration_seasons", "median"),
            max_duration=("duration_seasons", "max"),
        )
    )
    return summary


def enso_persistence_index(blocks: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula índices comparativos de persistência ENSO por série.

    Definições usadas no script:
        ENSO_event_fraction = fração dos trimestres classificados como EN ou LN.
        ENSO_event_count = número de blocos EN/LN.
        ENSO_mean_event_duration = duração média dos blocos EN/LN.
        ENSO_max_event_duration = duração máxima dos blocos EN/LN.
        ENSO_persistence_index = duração média dos eventos EN/LN ponderada
                                 pela fração de trimestres EN/LN.

    Esse índice composto é útil para comparar, no mesmo conjunto temporal,
    se uma série tende a produzir eventos ENSO mais persistentes/frequentes.
    """
    if blocks.empty:
        return pd.DataFrame()

    rows = []
    for (series, index_name), sub in blocks.groupby(["series", "index"]):
        total = sub["duration_seasons"].sum()
        enso = sub[sub["class"].isin(["EN", "LN"])]
        en = sub[sub["class"] == "EN"]
        ln = sub[sub["class"] == "LN"]

        enso_total = enso["duration_seasons"].sum() if not enso.empty else 0
        enso_count = len(enso)
        enso_mean = enso["duration_seasons"].mean() if enso_count > 0 else 0.0
        enso_median = enso["duration_seasons"].median() if enso_count > 0 else 0.0
        enso_max = enso["duration_seasons"].max() if enso_count > 0 else 0
        enso_fraction = enso_total / total if total > 0 else np.nan
        persistence_index = enso_mean * enso_fraction if total > 0 else np.nan

        rows.append({
            "series": series,
            "index": index_name,
            "total_seasons": int(total),
            "ENSO_event_seasons": int(enso_total),
            "ENSO_event_fraction": float(enso_fraction),
            "ENSO_event_count": int(enso_count),
            "ENSO_mean_event_duration": float(enso_mean),
            "ENSO_median_event_duration": float(enso_median),
            "ENSO_max_event_duration": int(enso_max),
            "ENSO_persistence_index": float(persistence_index),
            "EN_event_count": int(len(en)),
            "EN_mean_duration": float(en["duration_seasons"].mean()) if len(en) > 0 else 0.0,
            "LN_event_count": int(len(ln)),
            "LN_mean_duration": float(ln["duration_seasons"].mean()) if len(ln) > 0 else 0.0,
        })

    return pd.DataFrame(rows)


def comparative_persistence_index(persistence: pd.DataFrame, ref_label: str, cmp_label: str) -> pd.DataFrame:
    """Compara índices de persistência entre comparação e referência."""
    if persistence.empty:
        return pd.DataFrame()

    ref = persistence[persistence["series"] == ref_label]
    cmp_ = persistence[persistence["series"] == cmp_label]

    if ref.empty or cmp_.empty:
        return pd.DataFrame()

    ref = ref.iloc[0]
    cmp_ = cmp_.iloc[0]

    metrics = [
        "ENSO_event_fraction",
        "ENSO_event_count",
        "ENSO_mean_event_duration",
        "ENSO_median_event_duration",
        "ENSO_max_event_duration",
        "ENSO_persistence_index",
        "EN_event_count",
        "EN_mean_duration",
        "LN_event_count",
        "LN_mean_duration",
    ]

    rows = []
    for metric in metrics:
        ref_v = ref[metric]
        cmp_v = cmp_[metric]
        diff = cmp_v - ref_v
        ratio = cmp_v / ref_v if pd.notna(ref_v) and ref_v != 0 else np.nan
        rows.append({
            "metric": metric,
            "reference_value": ref_v,
            "comparison_value": cmp_v,
            "difference_comparison_minus_reference": diff,
            "ratio_comparison_over_reference": ratio,
        })

    return pd.DataFrame(rows)


# ============================================================
# 6. MATRIZ DE CONFUSÃO E MÉTRICAS DE SKILL
# ============================================================

def confusion_matrix(aligned: pd.DataFrame) -> pd.DataFrame:
    """Matriz de confusão EN/N/LN com linhas=referência e colunas=comparação."""
    matrix = pd.crosstab(aligned["class_ref"], aligned["class_cmp"])
    matrix = matrix.reindex(index=CLASS_ORDER, columns=CLASS_ORDER, fill_value=0)
    return matrix


def safe_div(num: float, den: float) -> float:
    """Divisão segura."""
    if den == 0 or pd.isna(den):
        return np.nan
    return num / den


def compute_skill_metrics(matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula métricas categóricas de skill a partir da matriz de confusão.

    Inclui:
        - Accuracy geral
        - Cohen's kappa / Heidke Skill Score multicategoria
        - Precision, recall, specificity, F1 e Peirce/TSS por classe
        - ETS por classe
        - médias macro das métricas por classe
    """
    mat = matrix.reindex(index=CLASS_ORDER, columns=CLASS_ORDER, fill_value=0).astype(float)
    arr = mat.values
    total = arr.sum()
    correct = np.trace(arr)

    rows = []

    accuracy = safe_div(correct, total)
    row_totals = arr.sum(axis=1)
    col_totals = arr.sum(axis=0)
    expected_correct = safe_div(np.sum(row_totals * col_totals), total)

    # Cohen's kappa / Heidke skill score for multicategory contingency table.
    po = accuracy
    pe = safe_div(expected_correct, total)
    kappa = safe_div(po - pe, 1 - pe)

    # Formula equivalent for multicategory Heidke Skill Score.
    hss_num = 2 * (total * correct - np.sum(row_totals * col_totals))
    hss_den = total**2 - np.sum(row_totals * col_totals)
    hss = safe_div(hss_num, hss_den)

    rows.append({"scope": "overall", "class": "ALL", "metric": "n", "value": total})
    rows.append({"scope": "overall", "class": "ALL", "metric": "accuracy", "value": accuracy})
    rows.append({"scope": "overall", "class": "ALL", "metric": "cohen_kappa", "value": kappa})
    rows.append({"scope": "overall", "class": "ALL", "metric": "heidke_skill_score", "value": hss})

    class_metric_rows = []

    for i, cls in enumerate(CLASS_ORDER):
        tp = arr[i, i]
        fn = row_totals[i] - tp
        fp = col_totals[i] - tp
        tn = total - tp - fn - fp

        precision = safe_div(tp, tp + fp)
        recall = safe_div(tp, tp + fn)  # POD
        specificity = safe_div(tn, tn + fp)
        false_alarm_ratio = safe_div(fp, tp + fp)
        false_alarm_rate = safe_div(fp, fp + tn)  # POFD
        f1 = safe_div(2 * precision * recall, precision + recall)
        peirce_tss = recall - false_alarm_rate if pd.notna(recall) and pd.notna(false_alarm_rate) else np.nan

        random_hits = safe_div((tp + fn) * (tp + fp), total)
        ets = safe_div(tp - random_hits, tp + fn + fp - random_hits)

        metrics = {
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "tn": tn,
            "precision": precision,
            "recall_pod": recall,
            "specificity": specificity,
            "false_alarm_ratio": false_alarm_ratio,
            "false_alarm_rate_pofd": false_alarm_rate,
            "f1": f1,
            "peirce_tss": peirce_tss,
            "equitable_threat_score": ets,
        }

        for metric, value in metrics.items():
            class_metric_rows.append({"scope": "per_class", "class": cls, "metric": metric, "value": value})

    rows.extend(class_metric_rows)

    df = pd.DataFrame(rows)

    # Macro averages for selected per-class metrics.
    macro_metrics = [
        "precision", "recall_pod", "specificity", "f1",
        "peirce_tss", "equitable_threat_score",
    ]
    for metric in macro_metrics:
        vals = df[(df["scope"] == "per_class") & (df["metric"] == metric)]["value"]
        df = pd.concat([
            df,
            pd.DataFrame([{
                "scope": "macro_average",
                "class": "ALL",
                "metric": f"macro_{metric}",
                "value": vals.mean(skipna=True),
            }])
        ], ignore_index=True)

    return df


# ============================================================
# 7. FIGURA 1 — SÉRIE TEMPORAL + DIVERGÊNCIAS
# ============================================================

def compute_y_limits_from_aligned(aligned: pd.DataFrame) -> Tuple[float, float, float]:
    vals = pd.concat([aligned["value_ref"], aligned["value_cmp"]], ignore_index=True).dropna()
    if vals.empty:
        ymin_data, ymax_data = -3.0, 3.0
    else:
        ymin_data = min(float(vals.min()), COLD_THRESHOLD)
        ymax_data = max(float(vals.max()), WARM_THRESHOLD)
    span = ymax_data - ymin_data
    if span <= 0:
        span = 1.0
    ymin_plot = ymin_data - 0.10 * span
    ymax_plot = ymax_data + TOP_BAND_FRACTION * span + 0.12 * span
    return ymin_plot, ymax_data, ymax_plot


def y_positions_for_divergences(ymin_plot: float, ymax_data: float) -> Dict[str, float]:
    span = ymax_data - ymin_plot
    if span <= 0:
        span = 1.0
    band_height = span * TOP_BAND_FRACTION
    step = band_height / (len(DIVERGENCE_ORDER) + 1)
    base = ymax_data + step
    return {code: base + i * step for i, code in enumerate(DIVERGENCE_ORDER)}


def create_time_series_figure(aligned: pd.DataFrame, ref_label: str, cmp_label: str) -> go.Figure:
    """Série temporal com pontos de divergência no topo."""
    ymin_plot, ymax_data, ymax_plot = compute_y_limits_from_aligned(aligned)
    y_positions = y_positions_for_divergences(ymin_plot, ymax_data)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=aligned["time"], y=aligned["value_ref"], mode="lines", name=ref_label,
        line=dict(color=REFERENCE_LINE_COLOR, width=REFERENCE_LINE_WIDTH),
        customdata=np.stack([aligned["year"], aligned["month"], aligned["season"], aligned["class_ref"]], axis=-1),
        hovertemplate=(
            "Data=%{x|%Y-%m}<br>Valor=%{y:.3f}<br>"
            "Ano=%{customdata[0]}<br>Mês=%{customdata[1]}<br>"
            "Trimestre=%{customdata[2]}<br>Classe=%{customdata[3]}"
            f"<extra>{ref_label}</extra>"
        ),
    ))

    fig.add_trace(go.Scatter(
        x=aligned["time"], y=aligned["value_cmp"], mode="lines", name=cmp_label,
        line=dict(color=COMPARISON_LINE_COLOR, width=COMPARISON_LINE_WIDTH),
        customdata=np.stack([aligned["year"], aligned["month"], aligned["season"], aligned["class_cmp"]], axis=-1),
        hovertemplate=(
            "Data=%{x|%Y-%m}<br>Valor=%{y:.3f}<br>"
            "Ano=%{customdata[0]}<br>Mês=%{customdata[1]}<br>"
            "Trimestre=%{customdata[2]}<br>Classe=%{customdata[3]}"
            f"<extra>{cmp_label}</extra>"
        ),
    ))

    xmin, xmax = aligned["time"].min(), aligned["time"].max()

    if SHOW_THRESHOLDS:
        fig.add_trace(go.Scatter(
            x=[xmin, xmax], y=[WARM_THRESHOLD, WARM_THRESHOLD], mode="lines",
            name=f"Limiar EN (+{WARM_THRESHOLD})",
            line=dict(color="red", width=1, dash="dash"), hoverinfo="skip",
        ))
        fig.add_trace(go.Scatter(
            x=[xmin, xmax], y=[COLD_THRESHOLD, COLD_THRESHOLD], mode="lines",
            name=f"Limiar LN ({COLD_THRESHOLD})",
            line=dict(color="blue", width=1, dash="dash"), hoverinfo="skip",
        ))

    if SHOW_ZERO_LINE:
        fig.add_trace(go.Scatter(
            x=[xmin, xmax], y=[0, 0], mode="lines", name="Zero",
            line=dict(color="gray", width=1, dash="dot"), hoverinfo="skip",
        ))

    # Faixa superior.
    fig.add_shape(
        type="rect", xref="x", yref="y",
        x0=xmin, x1=xmax, y0=ymax_data, y1=ymax_plot,
        fillcolor="rgba(245,245,245,0.75)", line=dict(width=0), layer="below",
    )
    fig.add_annotation(
        x=xmin, y=ymax_plot,
        text="Divergências de classificação",
        showarrow=False, xanchor="left", yanchor="top", font=dict(size=12, color="gray"),
    )

    for code in DIVERGENCE_ORDER:
        subset = aligned[aligned["divergence_type"] == code].copy()
        y_value = y_positions[code]
        if subset.empty:
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(size=DIVERGENCE_MARKER_SIZE, color=DIVERGENCE_COLORS[code]),
                name=DIVERGENCE_LABELS[code],
            ))
            continue

        fig.add_trace(go.Scatter(
            x=subset["time"], y=np.full(len(subset), y_value), mode="markers",
            marker=dict(size=DIVERGENCE_MARKER_SIZE, color=DIVERGENCE_COLORS[code], line=dict(color="black", width=0.7)),
            name=DIVERGENCE_LABELS[code],
            customdata=np.stack([
                subset["year"], subset["month"], subset["season"],
                subset["class_ref"], subset["class_cmp"], subset["value_ref"], subset["value_cmp"],
            ], axis=-1),
            hovertemplate=(
                "Divergência: " + DIVERGENCE_LABELS[code] + "<br>"
                "Data=%{x|%Y-%m}<br>Ano=%{customdata[0]}<br>Mês=%{customdata[1]}<br>"
                "Trimestre=%{customdata[2]}<br>Classe ref.=%{customdata[3]}<br>Classe comp.=%{customdata[4]}<br>"
                "Valor ref.=%{customdata[5]:.3f}<br>Valor comp.=%{customdata[6]:.3f}"
                "<extra></extra>"
            ),
        ))

    n_div = int(aligned["is_divergent"].sum())
    n_common = len(aligned)
    subtitle = (
        f"{ref_label} vs {cmp_label} | Divergências: {n_div}/{n_common} trimestres comuns"
        # f"{ref_label} vs {cmp_label}<br>"
        # f"EN: sequência ≥{MIN_CONSECUTIVE_SEASONS} > {WARM_THRESHOLD}; "
        # f"LN: sequência ≥{MIN_CONSECUTIVE_SEASONS} < {COLD_THRESHOLD}; "
        # f"Divergências: {n_div}/{n_common} trimestres comuns"
    )

    fig.update_layout(
        title=dict(text=f"Série temporal e divergências ENSO<br><sup>{subtitle}</sup>", x=0.5),
        template=PLOT_TEMPLATE,
        width=FIGURE_WIDTH,
        height=TIME_SERIES_HEIGHT,
        xaxis=dict(title="Tempo"),
        yaxis=dict(title=f"Valor ({REFERENCE_INDEX}/{COMPARISON_INDEX}) + faixa superior de divergências", range=[ymin_plot, ymax_plot]),
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        margin=dict(l=70, r=35, t=165, b=60),
    )
    return fig


# ============================================================
# 8. FIGURA 2 — MATRIZ DE CONFUSÃO
# ============================================================

def create_confusion_figure(matrix: pd.DataFrame, ref_label: str, cmp_label: str) -> go.Figure:
    z = matrix.values
    total = z.sum()
    percent = np.where(total > 0, z / total * 100.0, np.nan)
    text = np.empty_like(z, dtype=object)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            text[i, j] = f"{z[i, j]}<br>{percent[i, j]:.1f}%" if total > 0 else str(z[i, j])

    fig = go.Figure(data=go.Heatmap(
        z=z, x=matrix.columns.tolist(), y=matrix.index.tolist(), colorscale="Blues",
        text=text, texttemplate="%{text}",
        hovertemplate="Referência=%{y}<br>Comparação=%{x}<br>n=%{z}<extra></extra>",
        colorbar=dict(title="n"),
    ))
    fig.update_layout(
        title=dict(text=f"Matriz de confusão ENSO<br><sup>Linhas: {ref_label} | Colunas: {cmp_label}</sup>", x=0.5),
        xaxis=dict(title="Classe da comparação"),
        yaxis=dict(title="Classe da referência", autorange="reversed"),
        width=850,
        height=CONFUSION_HEIGHT,
        template=PLOT_TEMPLATE,
        margin=dict(l=80, r=40, t=100, b=70),
    )
    return fig


# ============================================================
# 9. FIGURA 3 — MÉTRICAS DE SKILL
# ============================================================

def create_skill_figure(skill: pd.DataFrame) -> go.Figure:
    """Tabela interativa de métricas de skill."""
    display = skill.copy()
    display["value_fmt"] = display["value"].map(lambda v: f"{v:.4f}" if pd.notna(v) else "NaN")

    fig = go.Figure(data=[go.Table(
        header=dict(
            values=["Escopo", "Classe", "Métrica", "Valor"],
            fill_color="#d9eaf7",
            align="left",
            font=dict(size=12),
        ),
        cells=dict(
            values=[display["scope"], display["class"], display["metric"], display["value_fmt"]],
            fill_color="#f7fbff",
            align="left",
            font=dict(size=11),
            height=24,
        ),
    )])
    fig.update_layout(
        title=dict(text="Métricas categóricas de skill EN/LN/N", x=0.5),
        width=FIGURE_WIDTH,
        height=SKILL_HEIGHT,
        template=PLOT_TEMPLATE,
        margin=dict(l=30, r=30, t=80, b=30),
    )
    return fig


# ============================================================
# 10. FIGURA 4 — PERSISTÊNCIA
# ============================================================

def create_persistence_figure(blocks: pd.DataFrame, summary: pd.DataFrame, persistence: pd.DataFrame, persistence_cmp: pd.DataFrame) -> go.Figure:
    if blocks.empty:
        fig = go.Figure()
        fig.update_layout(title="Sem blocos para análise de persistência")
        return fig

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(
            "Duração dos blocos EN/LN/N",
            "Resumo por classe",
            # "Índices de persistência ENSO",
            # "Comparativo: comparação - referência",
        ),
        # specs=[[{"type": "xy"}, {"type": "table"}], [{"type": "table"}, {"type": "table"}]],
        specs=[[{"type": "xy"}, {"type": "table"}]],
        column_widths=[0.58, 0.42],
        # row_heights=[0.58, 0.42],
        horizontal_spacing=0.08,
        vertical_spacing=0.16,
    )

    class_colors = {"EN": "red", "N": "gray", "LN": "blue"}
    plot_blocks = blocks.copy()
    plot_blocks["block_label"] = (
        plot_blocks["series"] + " | " + plot_blocks["class"] + " | " +
        plot_blocks["start_season"].astype(str) + "/" + plot_blocks["start_year"].astype(str) + "–" +
        plot_blocks["end_season"].astype(str) + "/" + plot_blocks["end_year"].astype(str)
    )

    for cls in CLASS_ORDER:
        sub = plot_blocks[plot_blocks["class"] == cls]
        if sub.empty:
            continue
        fig.add_trace(go.Bar(
            x=sub["block_label"], y=sub["duration_seasons"], name=f"Blocos {cls}",
            marker_color=class_colors.get(cls, "gray"),
            customdata=np.stack([sub["series"], sub["class"], sub["start_season"], sub["start_year"], sub["end_season"], sub["end_year"]], axis=-1),
            hovertemplate=(
                "Série=%{customdata[0]}<br>Classe=%{customdata[1]}<br>"
                "Início=%{customdata[2]}/%{customdata[3]}<br>Fim=%{customdata[4]}/%{customdata[5]}<br>"
                "Duração=%{y} trimestres<extra></extra>"
            ),
        ), row=1, col=1)

    add_table_trace(fig, summary, row=1, col=2, float_cols=["mean_duration", "median_duration"])
    # add_table_trace(fig, persistence, row=2, col=1, float_cols=["ENSO_event_fraction", "ENSO_mean_event_duration", "ENSO_median_event_duration", "ENSO_persistence_index", "EN_mean_duration", "LN_mean_duration"])
    # add_table_trace(fig, persistence_cmp, row=2, col=2, float_cols=["reference_value", "comparison_value", "difference_comparison_minus_reference", "ratio_comparison_over_reference"])

    fig.update_layout(
        title=dict(text="Persistência e duração dos regimes ENSO", x=0.5),
        template=PLOT_TEMPLATE,
        width=FIGURE_WIDTH,
        height=PERSISTENCE_HEIGHT,
        barmode="group",
        margin=dict(l=60, r=40, t=120, b=210),
    )
    fig.update_xaxes(tickangle=45, row=1, col=1)
    fig.update_yaxes(title_text="Duração (trimestres)", row=1, col=1)
    return fig


def add_table_trace(fig: go.Figure, df: pd.DataFrame, row: int, col: int, float_cols: Optional[List[str]] = None) -> None:
    """Adiciona uma tabela genérica a um subplot."""
    if float_cols is None:
        float_cols = []

    if df is None or df.empty:
        fig.add_trace(go.Table(
            header=dict(values=["Sem dados"], fill_color="#d9eaf7", align="left"),
            cells=dict(values=[[""]], fill_color="#f7fbff", align="left"),
        ), row=row, col=col)
        return

    display = df.copy()
    for fc in float_cols:
        if fc in display.columns:
            display[fc] = display[fc].map(lambda v: f"{v:.4f}" if pd.notna(v) else "NaN")

    fig.add_trace(go.Table(
        header=dict(values=list(display.columns), fill_color="#d9eaf7", align="left", font=dict(size=10)),
        cells=dict(values=[display[c].tolist() for c in display.columns], fill_color="#f7fbff", align="left", font=dict(size=9), height=22),
    ), row=row, col=col)


# ============================================================
# 11. HTML FINAL
# ============================================================

def make_html(fig_ts: go.Figure, fig_conf: go.Figure, fig_skill: go.Figure, fig_persist: go.Figure) -> str:
    style = """
    <style>
        body { font-family: Arial, Helvetica, sans-serif; margin: 24px; color: #222; }
        h1 { color: #1f4e79; margin-bottom: 6px; }
        h2 { color: #1f4e79; border-bottom: 1px solid #ddd; padding-bottom: 4px; margin-top: 30px; }
        .note {
            background: #f4f8fb;
            border-left: 4px solid #1f77b4;
            padding: 10px 14px;
            margin: 12px 0 20px 0;
            font-size: 14px;
        }
    </style>
    """

    ts_html = to_html(fig_ts, include_plotlyjs="cdn", full_html=False, config={"responsive": True, "displaylogo": False})
    conf_html = to_html(fig_conf, include_plotlyjs=False, full_html=False, config={"responsive": True, "displaylogo": False})
    skill_html = to_html(fig_skill, include_plotlyjs=False, full_html=False, config={"responsive": True, "displaylogo": False})
    persist_html = to_html(fig_persist, include_plotlyjs=False, full_html=False, config={"responsive": True, "displaylogo": False})

    html = f"""
    <html>
    <head>
        <meta charset="utf-8">
        <title>{HTML_TITLE}</title>
        {style}
    </head>
    <body>
        <h1>{HTML_TITLE}</h1>
        <div class="note">
            Este HTML compara duas curvas, que podem ser ONI e/ou RONI. A classificação EN/LN/N é calculada separadamente
            para cada curva usando sequências completas de pelo menos {MIN_CONSECUTIVE_SEASONS} trimestres. Pontos coloridos
            no topo da série temporal indicam divergências entre referência e comparação. A geração de CSVs é controlada pela
            variável <b>SAVE_CSV</b> no início do script.
        </div>
        <h2>1. Série temporal e divergências</h2>
        {ts_html}
        <h2>2. Matriz de confusão EN/LN/N</h2>
        {conf_html}
        <h2>3. Métricas de skill</h2>
        {skill_html}
        <h2>4. Persistência e duração ENSO</h2>
        {persist_html}
    </body>
    </html>
    """
    return html


# ============================================================
# 12. MAIN
# ============================================================

def main() -> None:
    validate_config()

    ref_label = make_label(REFERENCE_DIRECTORY, REFERENCE_INDEX, REFERENCE_LABEL)
    cmp_label = make_label(COMPARISON_DIRECTORY, COMPARISON_INDEX, COMPARISON_LABEL)

    ref_df = add_classification(read_index_file(REFERENCE_DIRECTORY, REFERENCE_INDEX))
    cmp_df = add_classification(read_index_file(COMPARISON_DIRECTORY, COMPARISON_INDEX))

    aligned = align_reference_and_comparison(ref_df, cmp_df)
    matrix = confusion_matrix(aligned)
    skill = compute_skill_metrics(matrix)

    ref_blocks = extract_class_blocks(ref_df, ref_label, REFERENCE_INDEX)
    cmp_blocks = extract_class_blocks(cmp_df, cmp_label, COMPARISON_INDEX)
    blocks = pd.concat([ref_blocks, cmp_blocks], ignore_index=True)

    summary = persistence_summary(blocks)
    persistence = enso_persistence_index(blocks)
    persistence_cmp = comparative_persistence_index(persistence, ref_label, cmp_label)

    if SAVE_CSV:
        aligned.to_csv(OUTPUT_ALIGNED_DIAGNOSTIC_CSV, index=False)
        matrix.to_csv(OUTPUT_CONFUSION_MATRIX_CSV)
        skill.to_csv(OUTPUT_SKILL_METRICS_CSV, index=False)
        blocks.to_csv(OUTPUT_EVENT_BLOCKS_CSV, index=False)
        summary.to_csv(OUTPUT_PERSISTENCE_SUMMARY_CSV, index=False)
        persistence_cmp.to_csv(OUTPUT_PERSISTENCE_COMPARISON_CSV, index=False)

    fig_ts = create_time_series_figure(aligned, ref_label, cmp_label)
    fig_conf = create_confusion_figure(matrix, ref_label, cmp_label)
    fig_skill = create_skill_figure(skill)
    fig_persist = create_persistence_figure(blocks, summary, persistence, persistence_cmp)

    html = make_html(fig_ts, fig_conf, fig_skill, fig_persist)
    Path(OUTPUT_HTML).write_text(html, encoding="utf-8")

    print("Processamento concluído.")
    print(f"HTML gerado: {OUTPUT_HTML}")
    print(f"SAVE_CSV: {SAVE_CSV}")
    if SAVE_CSV:
        print("CSVs gerados:")
        print(f"  - {OUTPUT_ALIGNED_DIAGNOSTIC_CSV}")
        print(f"  - {OUTPUT_CONFUSION_MATRIX_CSV}")
        print(f"  - {OUTPUT_SKILL_METRICS_CSV}")
        print(f"  - {OUTPUT_EVENT_BLOCKS_CSV}")
        print(f"  - {OUTPUT_PERSISTENCE_SUMMARY_CSV}")
        print(f"  - {OUTPUT_PERSISTENCE_COMPARISON_CSV}")
    print(f"Referência: {ref_label}")
    print(f"Comparação: {cmp_label}")
    print(f"Trimestres comuns: {len(aligned)}")
    print(f"Trimestres divergentes: {int(aligned['is_divergent'].sum())}")


if __name__ == "__main__":
    main()
