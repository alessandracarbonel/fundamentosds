#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_oni_roni_directories_plotly.py

Cria um HTML interativo em Plotly para comparar séries ONI e/ou RONI geradas
pelos scripts anteriores. Cada diretório comparado deve conter os arquivos:

    - ONI_timeseries.csv
    - RONI_timeseries.csv

Estrutura do HTML para cada índice:
    Linha 1: série temporal de todos os diretórios.
    Linha 2: gráficos de dispersão entre uma série de referência e cada série comparada.
    Linha 3: tabela com métricas estatísticas sob cada gráfico de dispersão.

Métricas calculadas:
    - RMSE
    - MAE
    - Bias
    - Correlação de Pearson (r)
    - R²
    - Slope
    - Intercept
    - Desvio-padrão do erro
    - n

Requisitos:
    pip install pandas numpy plotly

Observação importante:
    Um arquivo HTML estático não pode, por segurança do navegador, abrir diretórios
    locais arbitrariamente depois de criado. Por isso, a seleção dos diretórios é feita
    no momento da execução deste script Python. O HTML gerado é interativo para zoom,
    pan, legenda, hover, download do gráfico etc.
"""

# ============================================================
# 1. CONFIGURAÇÃO DO USUÁRIO
# ============================================================

# ------------------------------------------------------------
# Modo de escolha dos diretórios
# ------------------------------------------------------------
# True  -> abre janelas para seleção interativa de diretórios usando tkinter.
# False -> usa a lista DIRECTORIES definida abaixo.
USE_DIRECTORY_DIALOG = False

# Lista de diretórios a comparar, caso USE_DIRECTORY_DIALOG = False.
# Mínimo: 2 diretórios. Máximo: 4 diretórios.
# O primeiro diretório da lista será usado como referência, a menos que
# REFERENCE_DIR_INDEX seja alterado.
DIRECTORIES = [
    r"./resultados/ersstv5",
    r"./resultados/era5",
    r"./resultados/era5_56_25",
    # r"./teste_3",
]

# Índice do diretório de referência dentro da lista DIRECTORIES.
# Ex.: 0 significa que o primeiro diretório será a referência.
REFERENCE_DIR_INDEX = 0

# ------------------------------------------------------------
# Índices a processar
# ------------------------------------------------------------
# Opções:
#   ["ONI"]
#   ["RONI"]
#   ["ONI", "RONI"]
INDICES_TO_COMPARE = ["ONI", "RONI"]

# Arquivos esperados dentro de cada diretório.
INDEX_FILES = {
    "ONI": "ONI_timeseries.csv",
    "RONI": "RONI_timeseries.csv",
}

# Colunas esperadas em cada arquivo.
INDEX_VALUE_COLUMNS = {
    "ONI": "ONI",
    "RONI": "RONI",
}

# Colunas preferenciais para alinhar as séries.
# O script tenta year/month primeiro. Se não existirem, usa time.
MERGE_KEYS_PREFERENCE = ["year_month", "time"]

# ------------------------------------------------------------
# Saídas
# ------------------------------------------------------------
OUTPUT_HTML = "ONI_RONI_directory_comparison.html"
# OUTPUT_METRICS_CSV = "ONI_RONI_directory_comparison_metrics.csv"

# ------------------------------------------------------------
# Aparência do HTML/gráficos
# ------------------------------------------------------------
HTML_TITLE = "Comparação Interativa de ONI/RONI entre Diretórios"
FIGURE_HEIGHT_PER_INDEX = 900
FIGURE_WIDTH = 1400
PLOT_TEMPLATE = "plotly_white"

# Mostrar limiares ENSO +/-0.5 na série temporal.
SHOW_ENSO_THRESHOLDS = True

# Se True, remove linhas com NaN antes de comparar.
DROP_NA_BEFORE_COMPARISON = True

# Número de casas decimais nas tabelas de métricas.
METRICS_DECIMALS = 4

# Prefixo usado no nome das séries quando o diretório não tiver nome claro.
DEFAULT_LABEL_PREFIX = "Serie"


# ============================================================
# 2. IMPORTS
# ============================================================

import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import to_html


# ============================================================
# 3. SELEÇÃO E VALIDAÇÃO DOS DIRETÓRIOS
# ============================================================

def select_directories_with_dialog(min_dirs: int = 2, max_dirs: int = 4) -> List[str]:
    """
    Seleciona diretórios usando interface gráfica tkinter.

    A seleção ocorre em janelas sucessivas. O usuário deve selecionar ao menos
    min_dirs e no máximo max_dirs diretórios. Cancelar após o mínimo encerra a seleção.
    """
    import tkinter as tk
    from tkinter import filedialog, messagebox

    root = tk.Tk()
    root.withdraw()

    selected = []

    for i in range(max_dirs):
        msg = (
            f"Selecione o diretório {i + 1} de até {max_dirs}.\n\n"
            f"Mínimo necessário: {min_dirs}.\n"
            "Após selecionar o mínimo, cancele para finalizar."
        )
        messagebox.showinfo("Selecionar diretório", msg)
        folder = filedialog.askdirectory(title=f"Selecionar diretório {i + 1}")

        if not folder:
            if len(selected) >= min_dirs:
                break
            messagebox.showwarning(
                "Seleção insuficiente",
                f"É necessário selecionar pelo menos {min_dirs} diretórios.",
            )
            continue

        selected.append(folder)

    root.destroy()

    if len(selected) < min_dirs:
        raise ValueError(f"Selecione pelo menos {min_dirs} diretórios.")

    return selected


def normalize_directories() -> List[Path]:
    """
    Obtém e valida a lista de diretórios.
    """
    if USE_DIRECTORY_DIALOG:
        dirs = select_directories_with_dialog(min_dirs=2, max_dirs=4)
    else:
        dirs = DIRECTORIES

    if len(dirs) < 2:
        raise ValueError("É necessário informar pelo menos 2 diretórios.")

    if len(dirs) > 4:
        raise ValueError("É permitido comparar no máximo 4 diretórios.")

    paths = [Path(d).expanduser().resolve() for d in dirs]

    for p in paths:
        if not p.exists() or not p.is_dir():
            raise FileNotFoundError(f"Diretório não encontrado: {p}")

    if REFERENCE_DIR_INDEX < 0 or REFERENCE_DIR_INDEX >= len(paths):
        raise ValueError("REFERENCE_DIR_INDEX está fora do intervalo de diretórios.")

    return paths


def get_directory_label(path: Path, existing_labels: List[str]) -> str:
    """
    Cria rótulo curto para o diretório.
    """
    label = path.name.strip() or f"{DEFAULT_LABEL_PREFIX}_{len(existing_labels) + 1}"

    original_label = label
    counter = 2
    while label in existing_labels:
        label = f"{original_label}_{counter}"
        counter += 1

    return label


# ============================================================
# 4. LEITURA E PREPARAÇÃO DOS DADOS
# ============================================================

def read_index_file(directory: Path, index_name: str) -> pd.DataFrame:
    """
    Lê ONI_timeseries.csv ou RONI_timeseries.csv em um diretório.
    """
    if index_name not in INDEX_FILES:
        raise KeyError(f"Índice não configurado em INDEX_FILES: {index_name}")

    file_path = directory / INDEX_FILES[index_name]

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {file_path}")

    df = pd.read_csv(file_path)

    value_col = INDEX_VALUE_COLUMNS[index_name]

    if value_col not in df.columns:
        raise KeyError(
            f"Coluna '{value_col}' não encontrada em {file_path}. "
            f"Colunas disponíveis: {list(df.columns)}"
        )

    # Garantir tempo.
    if "time" in df.columns:
        df["time"] = pd.to_datetime(df["time"])
    elif {"year", "month"}.issubset(df.columns):
        df["time"] = pd.to_datetime(dict(year=df["year"], month=df["month"], day=15))
    else:
        raise KeyError(
            f"O arquivo {file_path} precisa conter 'time' ou as colunas 'year' e 'month'."
        )

    # Garantir year/month.
    if "year" not in df.columns:
        df["year"] = df["time"].dt.year
    if "month" not in df.columns:
        df["month"] = df["time"].dt.month

    if "season" not in df.columns:
        season_by_month = {
            1: "DJF", 2: "JFM", 3: "FMA", 4: "MAM", 5: "AMJ", 6: "MJJ",
            7: "JJA", 8: "JAS", 9: "ASO", 10: "SON", 11: "OND", 12: "NDJ",
        }
        df["season"] = df["month"].map(season_by_month)

    df = df[["time", "year", "month", "season", value_col]].copy()
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")
    df = df.sort_values(["year", "month"]).reset_index(drop=True)

    return df


def read_all_data(directories: List[Path], index_name: str) -> Dict[str, pd.DataFrame]:
    """
    Lê o arquivo do índice escolhido em todos os diretórios.
    """
    data = {}
    labels = []

    for directory in directories:
        label = get_directory_label(directory, labels)
        labels.append(label)
        data[label] = read_index_file(directory, index_name)

    return data


def select_reference_label(labels: List[str]) -> str:
    """
    Retorna o rótulo da série de referência.
    """
    return labels[REFERENCE_DIR_INDEX]


def align_pair(
    reference_df: pd.DataFrame,
    comparison_df: pd.DataFrame,
    index_name: str,
    ref_label: str,
    comp_label: str,
) -> pd.DataFrame:
    """
    Alinha duas séries por year/month/season.
    """
    value_col = INDEX_VALUE_COLUMNS[index_name]

    ref = reference_df.rename(columns={value_col: f"{index_name}_{ref_label}"})
    comp = comparison_df.rename(columns={value_col: f"{index_name}_{comp_label}"})

    merged = pd.merge(
        ref[["time", "year", "month", "season", f"{index_name}_{ref_label}" ]],
        comp[["year", "month", "season", f"{index_name}_{comp_label}" ]],
        on=["year", "month", "season"],
        how="inner",
    )

    if DROP_NA_BEFORE_COMPARISON:
        merged = merged.dropna(subset=[f"{index_name}_{ref_label}", f"{index_name}_{comp_label}"])

    return merged


# ============================================================
# 5. MÉTRICAS ESTATÍSTICAS
# ============================================================

def compute_pair_metrics(
    df: pd.DataFrame,
    ref_col: str,
    comp_col: str,
    index_name: str,
    ref_label: str,
    comp_label: str,
) -> Dict[str, object]:
    """
    Calcula as mesmas estatísticas usadas na validação anterior.
    """
    valid = df[[ref_col, comp_col]].dropna()
    n = len(valid)

    if n == 0:
        return {
            "index": index_name,
            "reference": ref_label,
            "comparison": comp_label,
            "n": 0,
            "rmse": np.nan,
            "mae": np.nan,
            "bias": np.nan,
            "r": np.nan,
            "r2": np.nan,
            "slope": np.nan,
            "intercept": np.nan,
            "error_std": np.nan,
        }

    ref = valid[ref_col].to_numpy(float)
    comp = valid[comp_col].to_numpy(float)
    diff = comp - ref

    rmse = np.sqrt(np.mean(diff ** 2))
    mae = np.mean(np.abs(diff))
    bias = np.mean(diff)
    error_std = np.std(diff, ddof=1) if n > 1 else np.nan

    if n > 1:
        r = np.corrcoef(ref, comp)[0, 1]
        r2 = r ** 2
        try:
            slope, intercept = np.polyfit(ref, comp, 1)
        except Exception:
            slope, intercept = np.nan, np.nan
    else:
        r = r2 = slope = intercept = np.nan

    return {
        "index": index_name,
        "reference": ref_label,
        "comparison": comp_label,
        "n": int(n),
        "rmse": float(rmse),
        "mae": float(mae),
        "bias": float(bias),
        "r": float(r),
        "r2": float(r2),
        "slope": float(slope),
        "intercept": float(intercept),
        "error_std": float(error_std),
    }


def format_metric(value, decimals: int = METRICS_DECIMALS) -> str:
    """
    Formata número para tabela Plotly.
    """
    if value is None or pd.isna(value):
        return "NaN"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return f"{float(value):.{decimals}f}"


# ============================================================
# 6. FIGURA PLOTLY
# ============================================================

def create_index_figure(index_name: str, data: Dict[str, pd.DataFrame]) -> Tuple[go.Figure, pd.DataFrame]:
    """
    Cria figura Plotly com 3 linhas para um índice.
    """
    labels = list(data.keys())
    ref_label = select_reference_label(labels)
    comp_labels = [label for label in labels if label != ref_label]
    n_comp = len(comp_labels)

    if n_comp < 1:
        raise ValueError("É necessário pelo menos um diretório para comparar contra a referência.")

    subplot_titles = [
        f"{index_name}: séries temporais",
        *[f"Dispersão: {label} vs {ref_label}" for label in comp_labels],
        *[f"Métricas: {label} vs {ref_label}" for label in comp_labels],
    ]

    specs_row1 = [{"type": "xy", "colspan": n_comp}] + [None] * (n_comp - 1)
    specs_row2 = [{"type": "xy"} for _ in range(n_comp)]
    specs_row3 = [{"type": "table"} for _ in range(n_comp)]

    fig = make_subplots(
        rows=3,
        cols=n_comp,
        specs=[specs_row1, specs_row2, specs_row3],
        row_heights=[0.45, 0.35, 0.20],
        vertical_spacing=0.10,
        horizontal_spacing=0.075,
        subplot_titles=subplot_titles,
    )

    value_col = INDEX_VALUE_COLUMNS[index_name]

    # Linha 1: séries temporais de todos os diretórios.
    for label, df in data.items():
        fig.add_trace(
            go.Scatter(
                x=df["time"],
                y=df[value_col],
                mode="lines",
                name=f"{index_name} - {label}",
                hovertemplate="%{x|%Y-%m}<br>valor=%{y:.3f}<extra>" + label + "</extra>",
            ),
            row=1,
            col=1,
        )

    if SHOW_ENSO_THRESHOLDS:
        # Linhas de referência no painel temporal.
        all_times = pd.concat([df["time"] for df in data.values()])
        xmin, xmax = all_times.min(), all_times.max()
        for yval, color, name in [(0.5, "red", "+0.5"), (-0.5, "blue", "-0.5")]:
            fig.add_trace(
                go.Scatter(
                    x=[xmin, xmax],
                    y=[yval, yval],
                    mode="lines",
                    line=dict(color=color, dash="dash", width=1),
                    name=name,
                    showlegend=True,
                    hoverinfo="skip",
                ),
                row=1,
                col=1,
            )

    metrics_records = []

    # Linhas 2 e 3: dispersão e métricas.
    reference_df = data[ref_label]

    for i, comp_label in enumerate(comp_labels, start=1):
        comp_df = data[comp_label]
        aligned = align_pair(reference_df, comp_df, index_name, ref_label, comp_label)

        ref_col = f"{index_name}_{ref_label}"
        comp_col = f"{index_name}_{comp_label}"

        metrics = compute_pair_metrics(aligned, ref_col, comp_col, index_name, ref_label, comp_label)
        metrics_records.append(metrics)

        # Scatter.
        fig.add_trace(
            go.Scatter(
                x=aligned[ref_col],
                y=aligned[comp_col],
                mode="markers",
                name=f"{comp_label} vs {ref_label}",
                marker=dict(size=7, opacity=0.75),
                customdata=np.stack(
                    [aligned["year"], aligned["month"], aligned["season"]],
                    axis=-1,
                ),
                hovertemplate=(
                    "Ano=%{customdata[0]}<br>"
                    "Mês=%{customdata[1]}<br>"
                    "Estação=%{customdata[2]}<br>"
                    f"{ref_label}=%{{x:.3f}}<br>"
                    f"{comp_label}=%{{y:.3f}}"
                    "<extra></extra>"
                ),
                showlegend=False,
            ),
            row=2,
            col=i,
        )

        # Linha 1:1 no scatter.
        if not aligned.empty:
            min_val = np.nanmin([aligned[ref_col].min(), aligned[comp_col].min()])
            max_val = np.nanmax([aligned[ref_col].max(), aligned[comp_col].max()])
            fig.add_trace(
                go.Scatter(
                    x=[min_val, max_val],
                    y=[min_val, max_val],
                    mode="lines",
                    line=dict(color="black", dash="dash", width=1),
                    showlegend=False,
                    hoverinfo="skip",
                ),
                row=2,
                col=i,
            )

        # Tabela de métricas.
        metric_names = [
            "n",
            "RMSE",
            "MAE",
            "Bias",
            "r",
            "R²",
            "Slope",
            "Intercept",
            "Erro std",
        ]
        metric_keys = [
            "n",
            "rmse",
            "mae",
            "bias",
            "r",
            "r2",
            "slope",
            "intercept",
            "error_std",
        ]
        metric_values = [format_metric(metrics[k]) for k in metric_keys]

        fig.add_trace(
            go.Table(
                header=dict(
                    values=["Métrica", "Valor"],
                    fill_color="#d9eaf7",
                    align="left",
                    font=dict(size=11),
                ),
                cells=dict(
                    values=[metric_names, metric_values],
                    fill_color="#f7fbff",
                    align="left",
                    font=dict(size=10),
                    height=22,
                ),
            ),
            row=3,
            col=i,
        )

        fig.update_xaxes(title_text=f"{ref_label}", row=2, col=i)
        fig.update_yaxes(title_text=f"{comp_label}", row=2, col=i)

    fig.update_xaxes(title_text="Tempo", row=1, col=1)
    fig.update_yaxes(title_text=f"{index_name} (°C)", row=1, col=1)

    fig.update_layout(
        title=dict(text=f"Comparação Interativa - {index_name}", x=0.5),
        template=PLOT_TEMPLATE,
        height=FIGURE_HEIGHT_PER_INDEX,
        width=FIGURE_WIDTH,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        hovermode="closest",
    )

    metrics_df = pd.DataFrame(metrics_records)
    return fig, metrics_df


# ============================================================
# 7. HTML FINAL
# ============================================================

def make_html_document(figures: Dict[str, go.Figure], metrics_df: pd.DataFrame) -> str:
    """
    Combina uma ou mais figuras Plotly em um único HTML.
    """
    parts = []

    style = """
    <style>
        body {
            font-family: Arial, Helvetica, sans-serif;
            margin: 24px;
            background: #ffffff;
            color: #222222;
        }
        h1 {
            color: #1f4e79;
            margin-bottom: 6px;
        }
        h2 {
            color: #1f4e79;
            border-bottom: 1px solid #dddddd;
            padding-bottom: 4px;
            margin-top: 30px;
        }
        .note {
            background: #f4f8fb;
            border-left: 4px solid #1f77b4;
            padding: 10px 14px;
            margin: 14px 0 22px 0;
            font-size: 14px;
        }
        .metrics-summary {
            font-size: 13px;
            margin-top: 20px;
        }
        table.summary-table {
            border-collapse: collapse;
            font-size: 12px;
        }
        table.summary-table th, table.summary-table td {
            border: 1px solid #dddddd;
            padding: 5px 8px;
        }
        table.summary-table th {
            background: #d9eaf7;
        }
    </style>
    """

    parts.append("<html><head><meta charset='utf-8'>")
    parts.append(f"<title>{HTML_TITLE}</title>")
    parts.append(style)
    parts.append("</head><body>")
    parts.append(f"<h1>{HTML_TITLE}</h1>")
    parts.append(
        "<div class='note'>"
        "Este HTML foi gerado a partir dos arquivos <b>ONI_timeseries.csv</b> e/ou "
        "<b>RONI_timeseries.csv</b> encontrados nos diretórios selecionados. "
        "A primeira série definida como referência é comparada contra as demais. "
        "Os gráficos são interativos: zoom, pan, hover, seleção via legenda e exportação."
        "</div>"
    )

    include_plotlyjs = True
    for index_name, fig in figures.items():
        parts.append(f"<h2>{index_name}</h2>")
        parts.append(
            to_html(
                fig,
                include_plotlyjs="cdn" if include_plotlyjs else False,
                full_html=False,
                config={"responsive": True, "displaylogo": False},
            )
        )
        include_plotlyjs = False

    if not metrics_df.empty:
        parts.append("<h2>Resumo consolidado das métricas</h2>")
        parts.append("<div class='metrics-summary'>")
        parts.append(metrics_df.to_html(index=False, classes="summary-table", float_format=lambda x: f"{x:.4f}"))
        parts.append("</div>")

    parts.append("</body></html>")
    return "\n".join(parts)


# ============================================================
# 8. MAIN
# ============================================================

def main() -> None:
    directories = normalize_directories()

    print("Diretórios selecionados:")
    for i, d in enumerate(directories):
        tag = " [REFERÊNCIA]" if i == REFERENCE_DIR_INDEX else ""
        print(f"  {i}: {d}{tag}")

    figures = {}
    all_metrics = []

    for index_name in INDICES_TO_COMPARE:
        if index_name not in INDEX_FILES:
            raise KeyError(f"Índice não configurado: {index_name}")

        print(f"\nProcessando índice: {index_name}")
        data = read_all_data(directories, index_name)
        fig, metrics_df = create_index_figure(index_name, data)
        figures[index_name] = fig
        all_metrics.append(metrics_df)

    metrics_all = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    # metrics_all.to_csv(OUTPUT_METRICS_CSV, index=False, float_format="%.6f")

    html = make_html_document(figures, metrics_all)
    Path(OUTPUT_HTML).write_text(html, encoding="utf-8")

    print("\nProcessamento concluído.")
    print(f"HTML gerado: {OUTPUT_HTML}")
    # print(f"Métricas consolidadas: {OUTPUT_METRICS_CSV}")


if __name__ == "__main__":
    main()
