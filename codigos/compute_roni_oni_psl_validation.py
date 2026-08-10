#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_roni_oni_psl_validation.py

Calcula RONI e ONI a partir de SST mensal em NetCDF, valida contra séries
oficiais NOAA/PSL em formato texto e gera gráficos comparativos.

Principais saídas:
  - RONI_timeseries.csv
  - ONI_timeseries.csv
  - RONI_validation_comparison.csv
  - ONI_validation_comparison.csv
  - validation_metrics.csv
  - RONI_validation_plot.png
  - ONI_validation_plot.png

Séries oficiais usadas na validação:
  - ONI : https://psl.noaa.gov/data/correlation/oni.data
  - RONI: https://psl.noaa.gov/data/timeseries/month/data/roni.data
"""

# ============================================================
# 1. CONFIGURAÇÕES DO USUÁRIO
# ============================================================

# Fonte da SST: "noaa" para ERSSTv5 ou "local" para ERA5/outro NetCDF
DATA_SOURCE = "local"
INPUT_FILE = "C:\\Users\\ALESSANDRA\\Desktop\\POS\\Ciencias_Dados\\Projeto_Indices_Climaticos\\calculo_indices_oni_roni\\dados_sst\\sst.mnmean.nc"

NOAA_URL = "https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
NOAA_FILE = "sst.mnmean.nc"

SST_VAR_NAME = "sst"
SST_UNITS = "auto"        # "auto", "kelvin" ou "celsius"

# RONI
RONI_CLIM_START = 1991
RONI_CLIM_END = 2020
LAT_MAX_TROP = 20
LAT_MIN_TROP = -20
RONI_SCALE_TO_NINO34_VARIANCE = True

# ONI CPC-style
ONI_FORCE_LAST_CLIM_START = 1996  # use None para automático

# Período exportado, aplicado somente ao salvar os CSVs
OUTPUT_START_YEAR = 1956
OUTPUT_END_YEAR = 2025
DROP_NA_INDEX_ROWS = False

# Validação PSL/NOAA
RUN_VALIDATION = True
OFFICIAL_ONI_PSL_URL = "https://psl.noaa.gov/data/correlation/oni.data"
OFFICIAL_RONI_PSL_URL = "https://psl.noaa.gov/data/timeseries/month/data/roni.data"
OFFICIAL_ONI_LOCAL_FILE = "official_oni_psl.data"
OFFICIAL_RONI_LOCAL_FILE = "official_roni_psl.data"
VALIDATION_FORCE_DOWNLOAD = True
VALIDATION_ROUND_CALCULATED_TO_1_DECIMAL = True
VALIDATION_SAVE_RAW_METRICS_TOO = True

# Figuras
MAKE_VALIDATION_PLOTS = True
PLOT_CALCULATED_VERSION = "raw"  # "raw" ou "rounded"
SAVE_PLOT_PNG = True
SAVE_PLOT_PDF = False
PLOT_DPI = 300
PLOT_FIGSIZE = (13, 5)
PLOT_SHOW_ZERO_LINE = True
PLOT_SHOW_ENSO_THRESHOLDS = True

# Arquivos de saída
RONI_OUTPUT = "RONI_timeseries.csv"
ONI_OUTPUT = "ONI_timeseries.csv"
RONI_VALIDATION_COMPARISON_OUTPUT = "RONI_validation_comparison.csv"
ONI_VALIDATION_COMPARISON_OUTPUT = "ONI_validation_comparison.csv"
VALIDATION_METRICS_OUTPUT = "validation_metrics.csv"
ONI_VALIDATION_PLOT_PNG = "ONI_validation_plot.png"
RONI_VALIDATION_PLOT_PNG = "RONI_validation_plot.png"
ONI_VALIDATION_PLOT_PDF = "ONI_validation_plot.pdf"
RONI_VALIDATION_PLOT_PDF = "RONI_validation_plot.pdf"

# ============================================================
# 2. IMPORTS
# ============================================================

import os
from typing import Dict, Iterable, List, Optional, Tuple
from urllib.request import urlretrieve

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# ============================================================
# 3. FUNÇÕES UTILITÁRIAS
# ============================================================

def download_file(url: str, output_file: str, force: bool = False) -> None:
    if os.path.exists(output_file) and not force:
        print(f"Arquivo já existe. Usando: {output_file}")
        return
    print(f"Baixando: {url}")
    urlretrieve(url, output_file)
    print(f"Arquivo salvo: {output_file}")


def infer_coord_name(obj: xr.Dataset, candidates: Iterable[str]) -> Optional[str]:
    for name in candidates:
        if name in obj.coords or name in obj.dims:
            return name
    return None


def standardize_sst_dataset(ds: xr.Dataset, varname: str) -> xr.DataArray:
    if varname not in ds.data_vars:
        raise KeyError(
            f"Variável '{varname}' não encontrada. Variáveis disponíveis: {list(ds.data_vars)}"
        )

    da = ds[varname]

    squeeze_dims = [d for d in da.dims if da.sizes[d] == 1]
    if squeeze_dims:
        print(f"Removendo dimensões de tamanho 1: {squeeze_dims}")
        da = da.squeeze(dim=squeeze_dims, drop=True)

    tmp = da.to_dataset(name=varname)
    time_name = infer_coord_name(tmp, ["time", "valid_time"])
    lat_name = infer_coord_name(tmp, ["lat", "latitude"])
    lon_name = infer_coord_name(tmp, ["lon", "longitude"])

    if time_name is None:
        raise ValueError("Não foi possível identificar a coordenada temporal: time/valid_time.")
    if lat_name is None:
        raise ValueError("Não foi possível identificar a coordenada de latitude: lat/latitude.")
    if lon_name is None:
        raise ValueError("Não foi possível identificar a coordenada de longitude: lon/longitude.")

    rename = {}
    if time_name != "time": rename[time_name] = "time"
    if lat_name != "lat": rename[lat_name] = "lat"
    if lon_name != "lon": rename[lon_name] = "lon"
    if rename:
        print(f"Renomeando coordenadas/dimensões: {rename}")
        da = da.rename(rename)

    da = da.transpose("time", "lat", "lon", ...)

    extra = [d for d in da.dims if d not in ["time", "lat", "lon"]]
    if extra:
        print(f"Dimensões extras removidas por média: {extra}")
        da = da.mean(dim=extra, skipna=True)

    da = da.assign_coords(lon=(((da["lon"] + 180) % 360) - 180))
    da = da.sortby("lon")
    da = da.sortby("lat")
    return da


def convert_units_to_celsius(da: xr.DataArray, units: str = "auto") -> xr.DataArray:
    units = units.lower()
    if units not in ["auto", "kelvin", "celsius"]:
        raise ValueError("SST_UNITS deve ser 'auto', 'kelvin' ou 'celsius'.")
    if units == "kelvin":
        out = da - 273.15
        out.attrs["units"] = "degC"
        print("Convertendo SST Kelvin -> Celsius.")
        return out
    if units == "celsius":
        da.attrs["units"] = "degC"
        print("Assumindo SST em Celsius.")
        return da

    sample = da.isel(time=slice(0, min(12, da.sizes["time"]))).mean(skipna=True)
    try:
        sample_value = float(sample.compute().values)
    except Exception:
        sample_value = float(sample.values)
    if sample_value > 100.0:
        print("Unidade detectada como Kelvin. Convertendo para Celsius.")
        out = da - 273.15
    else:
        print("Unidade detectada como Celsius.")
        out = da
    out.attrs["units"] = "degC"
    return out


def clean_sst_values(da: xr.DataArray) -> xr.DataArray:
    return da.where(np.isfinite(da)).where(np.abs(da) < 100.0)


def area_weighted_mean(da: xr.DataArray) -> xr.DataArray:
    weights = np.cos(np.deg2rad(da["lat"]))
    return da.weighted(weights).mean(dim=["lat", "lon"], skipna=True)


def monthly_climatology(da: xr.DataArray, start_year: int, end_year: int) -> xr.DataArray:
    sub = da.sel(time=slice(f"{start_year}-01-01", f"{end_year}-12-31"))
    if sub.time.size == 0:
        raise ValueError(f"Não há dados para climatologia {start_year}-{end_year}.")
    return sub.groupby("time.month").mean("time", skipna=True)


def add_cpc_season_labels(df: pd.DataFrame, time_col: str = "time") -> pd.DataFrame:
    season_by_month = {
        1: "DJF", 2: "JFM", 3: "FMA", 4: "MAM", 5: "AMJ", 6: "MJJ",
        7: "JJA", 8: "JAS", 9: "ASO", 10: "SON", 11: "OND", 12: "NDJ",
    }
    df = df.copy()
    df[time_col] = pd.to_datetime(df[time_col])
    df["year"] = df[time_col].dt.year
    df["month"] = df[time_col].dt.month
    df["season"] = df["month"].map(season_by_month)
    return df


def apply_output_period_filter(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["time"] = pd.to_datetime(df["time"])
    if OUTPUT_START_YEAR is not None:
        df = df[df["time"].dt.year >= OUTPUT_START_YEAR]
    if OUTPUT_END_YEAR is not None:
        df = df[df["time"].dt.year <= OUTPUT_END_YEAR]
    return df

# ============================================================
# 4. RONI
# ============================================================

def compute_roni(sst: xr.DataArray) -> xr.Dataset:
    print("\nCalculando RONI...")
    clim = monthly_climatology(sst, RONI_CLIM_START, RONI_CLIM_END)
    anom = sst.groupby("time.month") - clim

    nino34_anom = anom.sel(lat=slice(-5, 5), lon=slice(-170, -120))
    tropics_anom = anom.sel(lat=slice(LAT_MIN_TROP, LAT_MAX_TROP))

    nino34_ts = area_weighted_mean(nino34_anom)
    tropics_ts = area_weighted_mean(tropics_anom)
    rel_raw = nino34_ts - tropics_ts

    if RONI_SCALE_TO_NINO34_VARIANCE:
        ref_nino = nino34_ts.sel(time=slice(f"{RONI_CLIM_START}-01-01", f"{RONI_CLIM_END}-12-31"))
        ref_rel = rel_raw.sel(time=slice(f"{RONI_CLIM_START}-01-01", f"{RONI_CLIM_END}-12-31"))
        scale = ref_nino.std(skipna=True) / ref_rel.std(skipna=True)
        rel_scaled = rel_raw * scale
    else:
        scale = xr.DataArray(1.0)
        rel_scaled = rel_raw

    roni = rel_scaled.rolling(time=3, center=True).mean()
    roni.name = "RONI"

    out = xr.Dataset({
        "nino34_anomaly": nino34_ts,
        "tropical_mean_anomaly": tropics_ts,
        "relative_nino34_raw": rel_raw,
        "relative_nino34_scaled": rel_scaled,
        "RONI": roni,
    })
    out.attrs["roni_climatology"] = f"{RONI_CLIM_START}-{RONI_CLIM_END}"
    out.attrs["roni_scale_factor"] = float(scale.values)
    return out

# ============================================================
# 5. ONI CPC-style
# ============================================================

def first_aligned_clim_start(first_year: int) -> int:
    y = int(first_year)
    while y % 5 != 1:
        y += 1
    return y


def last_aligned_clim_start(last_year: int) -> int:
    y = int(last_year) - 29
    while y % 5 != 1:
        y -= 1
    return y


def oni_nominal_clim_start_for_year(year: int) -> int:
    block_start = ((int(year) - 1) // 5) * 5 + 1
    return block_start - 15


def build_oni_climatology_mapping(
    years: np.ndarray,
    data_first_year: int,
    data_last_year: int,
    force_last_clim_start: Optional[int] = None,
) -> Dict[int, Tuple[int, int]]:
    first_start = first_aligned_clim_start(data_first_year)
    last_start = last_aligned_clim_start(data_last_year)

    if force_last_clim_start is not None:
        forced_end = force_last_clim_start + 29
        if forced_end <= data_last_year:
            last_start = min(last_start, force_last_clim_start)
        else:
            print(
                f"Atenção: {force_last_clim_start}-{forced_end} não está completo; "
                "usando último período completo disponível."
            )

    if first_start + 29 > data_last_year:
        raise ValueError("A série não possui pelo menos 30 anos completos para ONI.")

    mapping = {}
    for y in years:
        nominal = oni_nominal_clim_start_for_year(int(y))
        cs = min(max(nominal, first_start), last_start)
        mapping[int(y)] = (int(cs), int(cs + 29))
    return mapping


def compute_oni(sst: xr.DataArray) -> xr.Dataset:
    print("\nCalculando ONI CPC-style...")
    nino34 = sst.sel(lat=slice(-5, 5), lon=slice(-170, -120))
    nino34_sst = area_weighted_mean(nino34)
    nino34_sst.name = "nino34_sst"

    years = np.unique(nino34_sst.time.dt.year.values)
    mapping = build_oni_climatology_mapping(
        years=years,
        data_first_year=int(years.min()),
        data_last_year=int(years.max()),
        force_last_clim_start=ONI_FORCE_LAST_CLIM_START,
    )

    periods = sorted(set(mapping.values()))
    print(f"Períodos climatológicos ONI usados: primeiro {periods[0][0]}-{periods[0][1]}, último {periods[-1][0]}-{periods[-1][1]}")

    clim_cache = {(cs, ce): monthly_climatology(nino34_sst, cs, ce) for cs, ce in periods}

    values, cs_list, ce_list, times = [], [], [], []
    for t in nino34_sst.time.values:
        ts = pd.Timestamp(t)
        cs, ce = mapping[int(ts.year)]
        clim = clim_cache[(cs, ce)]
        val = nino34_sst.sel(time=t) - clim.sel(month=int(ts.month))
        values.append(float(val.values))
        cs_list.append(cs)
        ce_list.append(ce)
        times.append(ts)

    time_index = pd.to_datetime(times)
    nino34_anom = xr.DataArray(values, coords={"time": time_index}, dims=["time"], name="nino34_anomaly")
    oni = nino34_anom.rolling(time=3, center=True).mean()
    oni.name = "ONI"

    return xr.Dataset({
        "nino34_sst": xr.DataArray(nino34_sst.values, coords={"time": time_index}, dims=["time"]),
        "nino34_anomaly": nino34_anom,
        "ONI": oni,
        "climatology_start": xr.DataArray(cs_list, coords={"time": time_index}, dims=["time"]),
        "climatology_end": xr.DataArray(ce_list, coords={"time": time_index}, dims=["time"]),
    })

# ============================================================
# 6. EXPORTAÇÃO
# ============================================================

def export_roni_csv(ds_roni: xr.Dataset, output_file: str) -> pd.DataFrame:
    df = ds_roni.to_dataframe().reset_index()
    df = add_cpc_season_labels(df, "time")
    df = apply_output_period_filter(df)
    if DROP_NA_INDEX_ROWS:
        df = df.dropna(subset=["RONI"])
    cols = ["time", "year", "month", "season", "nino34_anomaly", "tropical_mean_anomaly", "relative_nino34_raw", "relative_nino34_scaled", "RONI"]
    df = df[cols]
    df.to_csv(output_file, index=False, float_format="%.4f")
    print(f"Arquivo salvo: {output_file}")
    return df


def export_oni_csv(ds_oni: xr.Dataset, output_file: str) -> pd.DataFrame:
    df = ds_oni.to_dataframe().reset_index()
    df = add_cpc_season_labels(df, "time")
    df = apply_output_period_filter(df)
    if DROP_NA_INDEX_ROWS:
        df = df.dropna(subset=["ONI"])
    cols = ["time", "year", "month", "season", "nino34_sst", "nino34_anomaly", "ONI", "climatology_start", "climatology_end"]
    df = df[cols]
    df.to_csv(output_file, index=False, float_format="%.4f")
    print(f"Arquivo salvo: {output_file}")
    return df

# ============================================================
# 7. VALIDAÇÃO PSL/NOAA
# ============================================================

def _is_missing_psl_value(value: float) -> bool:
    return (not np.isfinite(value)) or value <= -90.0 or value >= 90.0


def read_psl_monthly_index_data(url: str, local_file: str, value_name: str, force_download: bool = False) -> pd.DataFrame:
    download_file(url, local_file, force=force_download)
    rows = []
    with open(local_file, "r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("/"):
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            try:
                year = int(float(parts[0]))
            except Exception:
                continue
            if year < 1800 or year > 2300:
                continue
            for month, token in enumerate(parts[1:13], start=1):
                token = token.replace("−", "-").replace("*", "")
                try:
                    value = float(token)
                except Exception:
                    value = np.nan
                if _is_missing_psl_value(value):
                    value = np.nan
                rows.append({"year": year, "month": month, value_name: value})

    if not rows:
        raise ValueError(f"Nenhuma linha mensal válida foi encontrada em {local_file}")

    df = pd.DataFrame(rows)
    df["time"] = pd.to_datetime(dict(year=df.year, month=df.month, day=15))
    df = add_cpc_season_labels(df, "time")
    df = df[["time", "year", "month", "season", value_name]].dropna(subset=[value_name])
    df = df.sort_values(["year", "month"]).reset_index(drop=True)
    print(f"Série oficial {value_name}: {df.year.min()}-{df.year.max()} ({len(df)} registros válidos).")
    return df


def compute_validation_metrics(comparison: pd.DataFrame, calc_col: str, official_col: str, index_name: str, mode: str) -> Dict[str, float]:
    df = comparison[[calc_col, official_col]].dropna()
    n = len(df)
    if n == 0:
        return {"index": index_name, "mode": mode, "n": 0, "rmse": np.nan, "mae": np.nan, "bias": np.nan, "r": np.nan, "r2": np.nan, "slope": np.nan, "intercept": np.nan, "error_std": np.nan}
    calc = df[calc_col].to_numpy(float)
    official = df[official_col].to_numpy(float)
    diff = calc - official
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    bias = float(np.mean(diff))
    error_std = float(np.std(diff, ddof=1)) if n > 1 else np.nan
    if n > 1:
        r = float(np.corrcoef(official, calc)[0, 1])
        r2 = r ** 2
        slope, intercept = np.polyfit(official, calc, 1)
        slope, intercept = float(slope), float(intercept)
    else:
        r = r2 = slope = intercept = np.nan
    return {"index": index_name, "mode": mode, "n": int(n), "rmse": rmse, "mae": mae, "bias": bias, "r": r, "r2": r2, "slope": slope, "intercept": intercept, "error_std": error_std}


def plot_validation_comparison(comparison: pd.DataFrame, index_name: str, calculated_col: str, official_col: str, output_png: str, output_pdf: str) -> None:
    if not MAKE_VALIDATION_PLOTS:
        return
    df = comparison.copy()
    df["time"] = pd.to_datetime(df["time"])
    raw_col = f"{calculated_col}_raw_for_validation"
    eval_col = f"{calculated_col}_eval"
    if PLOT_CALCULATED_VERSION.lower() == "rounded":
        calc_col = eval_col
        calc_label = f"{index_name} calculado (arred. 1 casa)"
    else:
        calc_col = raw_col
        calc_label = f"{index_name} calculado"
    plot_df = df[["time", calc_col, official_col]].dropna()
    if plot_df.empty:
        print(f"Sem dados para plotar {index_name}.")
        return

    fig, ax = plt.subplots(figsize=PLOT_FIGSIZE)
    ax.plot(plot_df.time, plot_df[official_col], color="black", linewidth=1.7, label=f"{index_name} oficial PSL/NOAA")
    ax.plot(plot_df.time, plot_df[calc_col], color="tab:red", linewidth=1.2, alpha=0.85, label=calc_label)
    if PLOT_SHOW_ZERO_LINE:
        ax.axhline(0, color="0.35", linewidth=0.8)
    if PLOT_SHOW_ENSO_THRESHOLDS:
        ax.axhline(0.5, color="tab:red", linestyle="--", linewidth=0.9, alpha=0.65, label="+0.5 °C")
        ax.axhline(-0.5, color="tab:blue", linestyle="--", linewidth=0.9, alpha=0.65, label="-0.5 °C")
    ax.set_title(f"Comparação entre {index_name} calculado e oficial PSL/NOAA ({plot_df.time.dt.year.min()}–{plot_df.time.dt.year.max()})", fontsize=13, fontweight="bold")
    ax.set_xlabel("Ano")
    ax.set_ylabel(f"{index_name} (°C)")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    if SAVE_PLOT_PNG:
        fig.savefig(output_png, dpi=PLOT_DPI, bbox_inches="tight")
        print(f"Figura salva: {output_png}")
    if SAVE_PLOT_PDF:
        fig.savefig(output_pdf, dpi=PLOT_DPI, bbox_inches="tight")
        print(f"Figura salva: {output_pdf}")
    plt.close(fig)


def validate_one_index(calculated_df: pd.DataFrame, official_df: pd.DataFrame, index_name: str, calculated_col: str, official_col: str, comparison_output: str) -> Tuple[List[Dict[str, float]], pd.DataFrame]:
    calc = calculated_df.copy()
    off = official_df.copy()
    calc["year"] = calc.year.astype(int)
    calc["month"] = calc.month.astype(int)
    off["year"] = off.year.astype(int)
    off["month"] = off.month.astype(int)

    comparison = pd.merge(calc, off[["year", "month", "season", official_col]], on=["year", "month", "season"], how="inner")
    comparison[f"{calculated_col}_raw_for_validation"] = comparison[calculated_col]
    if VALIDATION_ROUND_CALCULATED_TO_1_DECIMAL:
        comparison[f"{calculated_col}_eval"] = comparison[calculated_col].round(1)
        eval_col = f"{calculated_col}_eval"
        mode = "calculated_rounded_1_decimal"
    else:
        comparison[f"{calculated_col}_eval"] = comparison[calculated_col]
        eval_col = f"{calculated_col}_eval"
        mode = "calculated_raw"

    raw_col = f"{calculated_col}_raw_for_validation"
    comparison["difference_eval_minus_official"] = comparison[eval_col] - comparison[official_col]
    comparison["difference_raw_minus_official"] = comparison[raw_col] - comparison[official_col]

    output_cols = ["time", "year", "month", "season", raw_col, eval_col, official_col, "difference_eval_minus_official", "difference_raw_minus_official"]
    comparison[output_cols].to_csv(comparison_output, index=False, float_format="%.4f")
    print(f"Arquivo de comparação salvo: {comparison_output}")

    metrics = [compute_validation_metrics(comparison, eval_col, official_col, index_name, mode)]
    if VALIDATION_SAVE_RAW_METRICS_TOO:
        metrics.append(compute_validation_metrics(comparison, raw_col, official_col, index_name, "calculated_raw"))
    return metrics, comparison[output_cols].copy()


def run_validation(roni_df: pd.DataFrame, oni_df: pd.DataFrame) -> pd.DataFrame:
    print("\nIniciando validação contra séries oficiais PSL/NOAA...")
    official_roni = read_psl_monthly_index_data(OFFICIAL_RONI_PSL_URL, OFFICIAL_RONI_LOCAL_FILE, "RONI_official", VALIDATION_FORCE_DOWNLOAD)
    official_oni = read_psl_monthly_index_data(OFFICIAL_ONI_PSL_URL, OFFICIAL_ONI_LOCAL_FILE, "ONI_official", VALIDATION_FORCE_DOWNLOAD)

    metrics_all = []
    roni_metrics, roni_comp = validate_one_index(roni_df, official_roni, "RONI", "RONI", "RONI_official", RONI_VALIDATION_COMPARISON_OUTPUT)
    metrics_all.extend(roni_metrics)
    plot_validation_comparison(roni_comp, "RONI", "RONI", "RONI_official", RONI_VALIDATION_PLOT_PNG, RONI_VALIDATION_PLOT_PDF)

    oni_metrics, oni_comp = validate_one_index(oni_df, official_oni, "ONI", "ONI", "ONI_official", ONI_VALIDATION_COMPARISON_OUTPUT)
    metrics_all.extend(oni_metrics)
    plot_validation_comparison(oni_comp, "ONI", "ONI", "ONI_official", ONI_VALIDATION_PLOT_PNG, ONI_VALIDATION_PLOT_PDF)

    metrics_df = pd.DataFrame(metrics_all)
    metrics_df.to_csv(VALIDATION_METRICS_OUTPUT, index=False, float_format="%.6f")
    print(f"Arquivo de métricas salvo: {VALIDATION_METRICS_OUTPUT}")
    print(metrics_df)
    return metrics_df

# ============================================================
# 8. MAIN
# ============================================================

def load_input_dataset() -> xr.Dataset:
    if DATA_SOURCE.lower() == "noaa":
        download_file(NOAA_URL, NOAA_FILE, force=False)
        input_file = NOAA_FILE
    elif DATA_SOURCE.lower() == "local":
        input_file = INPUT_FILE
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Arquivo local não encontrado: {input_file}")
    else:
        raise ValueError("DATA_SOURCE deve ser 'noaa' ou 'local'.")
    print(f"\nAbrindo arquivo: {input_file}")
    return xr.open_dataset(input_file, decode_times=True)


def main() -> None:
    ds = load_input_dataset()
    sst = standardize_sst_dataset(ds, SST_VAR_NAME)
    sst = convert_units_to_celsius(sst, SST_UNITS)
    sst = clean_sst_values(sst)
    print("\nResumo da SST padronizada:")
    print(sst)

    ds_roni = compute_roni(sst)
    ds_oni = compute_oni(sst)

    roni_df = export_roni_csv(ds_roni, RONI_OUTPUT)
    oni_df = export_oni_csv(ds_oni, ONI_OUTPUT)

    if RUN_VALIDATION:
        run_validation(roni_df, oni_df)

    print("\nProcessamento finalizado com sucesso.")
    print("Arquivos principais:")
    print(f"  - {RONI_OUTPUT}")
    print(f"  - {ONI_OUTPUT}")
    if RUN_VALIDATION:
        print("Arquivos de validação:")
        print(f"  - {RONI_VALIDATION_COMPARISON_OUTPUT}")
        print(f"  - {ONI_VALIDATION_COMPARISON_OUTPUT}")
        print(f"  - {VALIDATION_METRICS_OUTPUT}")
        if MAKE_VALIDATION_PLOTS and SAVE_PLOT_PNG:
            print(f"  - {RONI_VALIDATION_PLOT_PNG}")
            print(f"  - {ONI_VALIDATION_PLOT_PNG}")


if __name__ == "__main__":
    main()
