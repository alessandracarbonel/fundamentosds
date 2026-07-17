# Cálculo e Avaliação dos Índices ONI e RONI (Monitoramento do ENOS)

Ferramenta de código aberto, escrita em **Python**, para o **cálculo**, a **validação** e a
**análise comparativa/interativa** dos índices **ONI** (*Oceanic Niño Index*) e
**RONI** (*Relative Oceanic Niño Index*), usados no monitoramento operacional do
El Niño–Oscilação Sul (ENOS).

O projeto foi concebido com foco em **rastreabilidade** e **reprodutibilidade**, alinhado
aos princípios **FAIR** (*Findable, Accessible, Interoperable, Reusable*) e aos "4R's" da
ciência de dados (reprodutibilidade, replicabilidade, reutilização e rastreabilidade/robustez).
Toda a configuração fica no **início de cada script** (*config-driven workflow*), o que
favorece experimentos controlados de sensibilidade metodológica.

---

## 📁 Estrutura do repositório

```
calculo_indices_oni_roni/
├── compute_roni_oni_psl_validation.py            # Cálculo dos índices + validação contra NOAA/PSL
├── compare_oni_roni_directories_plotly.py        # Comparação interativa entre diretórios/bases (Plotly)
├── generate_enso_cross_index_skill_persistence_html.py  # Skill, matriz de confusão e persistência ENSO
│
├── ONI_timeseries.csv                            # Saída gerada (série ONI)
├── RONI_timeseries.csv                           # Saída gerada (série RONI)
├── ONI_RONI_directory_comparison.html            # Saída interativa (comparação entre diretórios)
├── ENSO_cross_index_skill_persistence.html       # Saída interativa (skill/persistência)
│
├── oni.data                                      # (opcional) referência oficial NOAA/PSL — validação ONI
└── roni.data                                     # (opcional) referência oficial NOAA/PSL — validação RONI
```

---

## 🧩 Descrição dos scripts

### 1. `compute_roni_oni_psl_validation.py`
Núcleo do projeto. Realiza:
- **Pré-processamento** da TSM (temperatura da superfície do mar): padronização de coordenadas,
  conversão de longitude (0–360° → −180–180°), ordenação de latitude, conversão Kelvin → °C,
  limpeza de valores não finitos/implausíveis e **média espacial ponderada pelo cosseno da latitude**.
- **Cálculo do RONI**: anomalia da região Niño 3.4 menos a anomalia média tropical (20°S–20°N),
  com ajuste de variância e média móvel de 3 meses.
- **Cálculo do ONI (estilo CPC/NOAA)**: climatologia móvel de 30 anos, atualizada a cada 5 anos e
  centrada nos blocos quinquenais.
- **Validação automática** das séries contra as referências oficiais da NOAA/PSL.
- **Saídas**: `ONI_timeseries.csv`, `RONI_timeseries.csv` e gráficos comparativos (calculado × oficial).

### 2. `compare_oni_roni_directories_plotly.py`
Gera o HTML interativo `ONI_RONI_directory_comparison.html` a partir dos arquivos
`ONI_timeseries.csv` e/ou `RONI_timeseries.csv` encontrados nos **diretórios selecionados**
(ex.: uma pasta com resultados ERSST v5 e outra com resultados ERA5). A primeira série definida
é usada como **referência** e comparada contra as demais. Produz gráficos interativos (zoom, pan,
hover, seleção via legenda, exportação) e um **resumo consolidado de métricas**
(RMSE, MAE, viés, r, R², slope, intercept, desvio-padrão do erro).

### 3. `generate_enso_cross_index_skill_persistence_html.py`
Gera o HTML interativo `ENSO_cross_index_skill_persistence.html`, comparando **duas curvas**
(qualquer combinação ONI×ONI, ONI×RONI, RONI×ONI, RONI×RONI). Inclui:
- **Série temporal** com marcação, no topo, dos tipos de divergência de classificação entre
  referência e comparação;
- **Classificação EN/LN/N** por sequências completas de ≥ 5 trimestres consecutivos;
- **Matriz de confusão** El Niño / La Niña / Neutro;
- **Métricas de skill**: Heidke Skill Score (HSS), Peirce/True Skill Statistic (TSS) e
  Equitable Threat Score (ETS);
- **Índice de persistência** e duração dos regimes ENOS.

---

## ⚙️ Requisitos

- **Python** ≥ 3.9
- Bibliotecas:

```bash
pip install numpy pandas xarray netCDF4 matplotlib plotly
```

| Biblioteca | Uso principal | Scripts |
|---|---|---|
| `numpy`, `pandas` | manipulação numérica e de séries temporais | todos |
| `xarray` + `netCDF4` | leitura de dados TSM em NetCDF (ERSST v5 / ERA5) | `compute_roni_oni_psl_validation.py` |
| `matplotlib` | gráficos de validação (calculado × oficial) | `compute_roni_oni_psl_validation.py` |
| `plotly` | gráficos interativos (HTML) | `compare_*` e `generate_enso_*` |
| `tkinter` | seleção interativa de diretórios (opcional, `USE_DIRECTORY_DIALOG=True`) | `compare_oni_roni_directories_plotly.py` |

> O download automático do ERSST v5 e das séries oficiais NOAA/PSL usa `urllib` (biblioteca
> padrão do Python), portanto **não** requer `requests`.
>
> Dependências mínimas por script:
> - `compute_roni_oni_psl_validation.py` → `numpy pandas xarray netCDF4 matplotlib`
> - `compare_oni_roni_directories_plotly.py` → `pandas numpy plotly`
> - `generate_enso_cross_index_skill_persistence_html.py` → `pandas numpy plotly`

> Recomenda-se o uso de um ambiente isolado (`venv`, `conda` ou `uv`) para garantir
> reprodutibilidade das versões.

Exemplo com `conda`:
```bash
conda create -n oni_roni python=3.10 numpy pandas xarray netCDF4 matplotlib plotly -c conda-forge
conda activate oni_roni
```

---

## 🌊 Dados de entrada

- **Fonte primária:** ERSST v5 (Extended Reconstructed SST v5), ~2°×2°, global, mensal — com
  **download automático** a partir dos servidores da NOAA/PSL (`sst.mnmean.nc`).
- **Fonte alternativa:** arquivo local NetCDF (ex.: reanálise **ERA5**).
- **Referências de validação (opcionais):** arquivos oficiais em formato PSL
  (`oni.data` e `roni.data`) da NOAA/CPC.

---

## 🔧 Configurações de usuário

Todas as configurações ficam em um **bloco no início de cada script** (seção
`1. CONFIGURAÇÕES DO USUÁRIO`). Abaixo, os parâmetros exatos de cada script, com seus
valores padrão.

### `compute_roni_oni_psl_validation.py`

**Fonte e leitura da TSM**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `DATA_SOURCE` | Fonte da TSM: `"noaa"` (download automático do ERSST v5) ou `"local"` (arquivo NetCDF, ex.: ERA5) | `"local"` |
| `INPUT_FILE` | Caminho do arquivo NetCDF local (usado quando `DATA_SOURCE="local"`) | `...\dados_sst\sst.mnmean.nc` |
| `NOAA_URL` | URL do ERSST v5 na NOAA/PSL (usado quando `DATA_SOURCE="noaa"`) | `https://downloads.psl.noaa.gov/.../sst.mnmean.nc` |
| `NOAA_FILE` | Nome do arquivo local salvo no download automático | `"sst.mnmean.nc"` |
| `SST_VAR_NAME` | Nome da variável de TSM dentro do NetCDF | `"sst"` |
| `SST_UNITS` | Tratamento de unidades: `"auto"`, `"kelvin"` ou `"celsius"` | `"auto"` |

**RONI**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `RONI_CLIM_START` | Ano inicial da climatologia fixa do RONI | `1991` |
| `RONI_CLIM_END` | Ano final da climatologia fixa do RONI | `2020` |
| `LAT_MIN_TROP` / `LAT_MAX_TROP` | Faixa tropical de remoção de tendência do RONI (ajustável) | `-20` / `20` |
| `RONI_SCALE_TO_NINO34_VARIANCE` | Aplica o ajuste de variância à Niño 3.4 | `True` |

**ONI (estilo CPC)**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `ONI_FORCE_LAST_CLIM_START` | Início do último período climatológico de 30 anos (use `None` para automático) | `1996` (→ 1996–2025) |

**Período exportado (aplicado apenas ao salvar os CSVs)**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `OUTPUT_START_YEAR` | Ano inicial de exportação | `1956` |
| `OUTPUT_END_YEAR` | Ano final de exportação | `2025` |
| `DROP_NA_INDEX_ROWS` | Remove linhas com índice `NaN` na exportação | `False` |

**Validação contra NOAA/PSL**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `RUN_VALIDATION` | Ativa a validação contra as séries oficiais | `True` |
| `OFFICIAL_ONI_PSL_URL` | URL da série oficial ONI (PSL) | `https://psl.noaa.gov/data/correlation/oni.data` |
| `OFFICIAL_RONI_PSL_URL` | URL da série oficial RONI (PSL) | `https://psl.noaa.gov/data/timeseries/month/data/roni.data` |
| `OFFICIAL_ONI_LOCAL_FILE` / `OFFICIAL_RONI_LOCAL_FILE` | Nomes locais dos arquivos oficiais baixados | `official_oni_psl.data` / `official_roni_psl.data` |
| `VALIDATION_FORCE_DOWNLOAD` | Força novo download das séries oficiais | `True` |
| `VALIDATION_ROUND_CALCULATED_TO_1_DECIMAL` | Arredonda a série calculada para 1 casa na validação | `True` |
| `VALIDATION_SAVE_RAW_METRICS_TOO` | Salva também as métricas sem arredondamento | `True` |

**Figuras de validação**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `MAKE_VALIDATION_PLOTS` | Gera gráficos de validação | `True` |
| `PLOT_CALCULATED_VERSION` | Versão plotada: `"raw"` ou `"rounded"` | `"raw"` |
| `SAVE_PLOT_PNG` / `SAVE_PLOT_PDF` | Formatos de saída dos gráficos | `True` / `False` |
| `PLOT_DPI` | Resolução das figuras | `300` |
| `PLOT_FIGSIZE` | Tamanho da figura (polegadas) | `(13, 5)` |
| `PLOT_SHOW_ZERO_LINE` / `PLOT_SHOW_ENSO_THRESHOLDS` | Mostra linha zero e limiares ±0,5 °C | `True` / `True` |

> **Regiões fixas do método:** Niño 3.4 = 5°N–5°S, 170°W–120°W (índices calculados);
> faixa tropical de referência = 20°N–20°S (remoção de tendência do RONI).
>
> **Saídas do script:** `RONI_timeseries.csv`, `ONI_timeseries.csv`,
> `RONI_validation_comparison.csv`, `ONI_validation_comparison.csv`, `validation_metrics.csv`,
> `RONI_validation_plot.png`, `ONI_validation_plot.png` (nomes configuráveis nas variáveis
> `*_OUTPUT` / `*_PLOT_*`).

### `compare_oni_roni_directories_plotly.py`

**Seleção dos diretórios**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `USE_DIRECTORY_DIALOG` | `True` → janela interativa de seleção (tkinter); `False` → usa a lista `DIRECTORIES` | `False` |
| `DIRECTORIES` | Lista de diretórios a comparar (mín. 2, máx. 4), cada um contendo `ONI_timeseries.csv` e/ou `RONI_timeseries.csv` | `["./resultados/ersstv5", "./resultados/era5", "./resultados/era5_56_25"]` |
| `REFERENCE_DIR_INDEX` | Índice do diretório de referência dentro de `DIRECTORIES` (0 = primeiro) | `0` |

**Índices e leitura**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `INDICES_TO_COMPARE` | Índices a processar: `["ONI"]`, `["RONI"]` ou `["ONI","RONI"]` | `["ONI", "RONI"]` |
| `INDEX_FILES` | Mapeamento índice → nome do arquivo CSV | `{"ONI":"ONI_timeseries.csv", "RONI":"RONI_timeseries.csv"}` |
| `INDEX_VALUE_COLUMNS` | Mapeamento índice → coluna de valor no CSV | `{"ONI":"ONI", "RONI":"RONI"}` |
| `MERGE_KEYS_PREFERENCE` | Chaves preferenciais de alinhamento das séries | `["year_month", "time"]` |

**Saída e aparência**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `OUTPUT_HTML` | Nome do HTML interativo gerado | `"ONI_RONI_directory_comparison.html"` |
| `HTML_TITLE` | Título exibido no HTML | `"Comparação Interativa de ONI/RONI entre Diretórios"` |
| `FIGURE_HEIGHT_PER_INDEX` / `FIGURE_WIDTH` | Dimensões da figura (px) | `900` / `1400` |
| `PLOT_TEMPLATE` | Tema visual do Plotly | `"plotly_white"` |
| `SHOW_ENSO_THRESHOLDS` | Mostra os limiares ENSO ±0,5 na série temporal | `True` |
| `DROP_NA_BEFORE_COMPARISON` | Remove linhas com `NaN` antes de comparar | `True` |
| `METRICS_DECIMALS` | Casas decimais nas tabelas de métricas | `4` |
| `DEFAULT_LABEL_PREFIX` | Prefixo de rótulo quando o diretório não tem nome claro | `"Serie"` |

> **Observação:** um HTML estático não pode acessar diretórios locais por segurança do navegador;
> por isso a seleção dos diretórios ocorre **na execução** do script. O HTML gerado permanece
> interativo (zoom, pan, hover, legenda, exportação).

### `generate_enso_cross_index_skill_persistence_html.py`

**Séries de entrada**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `REFERENCE_DIRECTORY` | Diretório da série de referência | `"./resultados/ersstv5"` |
| `COMPARISON_DIRECTORY` | Diretório da série comparada | `"./resultados/ersstv5"` |
| `REFERENCE_INDEX` | Índice da referência: `"ONI"` ou `"RONI"` | `"ONI"` |
| `COMPARISON_INDEX` | Índice da comparação: `"ONI"` ou `"RONI"` | `"RONI"` |
| `REFERENCE_LABEL` / `COMPARISON_LABEL` | Rótulos opcionais (`None` → `"ÍNDICE (pasta)"`) | `None` / `None` |

> Permite qualquer combinação: ONI×ONI, ONI×RONI, RONI×ONI, RONI×RONI.

**Critérios ENSO**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `WARM_THRESHOLD` | Limiar para El Niño (°C) | `0.5` |
| `COLD_THRESHOLD` | Limiar para La Niña (°C) | `-0.5` |
| `MIN_CONSECUTIVE_SEASONS` | Nº mínimo de trimestres consecutivos para caracterizar evento | `5` |

**Saídas**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `OUTPUT_HTML` | Nome do HTML interativo gerado | `"ENSO_cross_index_skill_persistence.html"` |
| `SAVE_CSV` | `True` → gera HTML + CSVs diagnósticos; `False` → apenas HTML | `False` |

> Quando `SAVE_CSV=True`, são gerados também: `ENSO_cross_index_aligned_diagnostic.csv`,
> `..._confusion_matrix.csv`, `..._skill_metrics.csv`, `..._event_blocks.csv`,
> `..._persistence_summary.csv` e `..._persistence_comparison.csv`.

**Aparência (opcional)**

| Parâmetro | Descrição | Valor padrão |
|---|---|---|
| `HTML_TITLE` | Título do HTML | `"Diagnóstico ENSO: comparação entre ONI/RONI, skill e persistência"` |
| `PLOT_TEMPLATE` | Tema visual do Plotly | `"plotly_white"` |
| `FIGURE_WIDTH` | Largura das figuras (px) | `1550` |
| `TIME_SERIES_HEIGHT` / `CONFUSION_HEIGHT` / `SKILL_HEIGHT` / `PERSISTENCE_HEIGHT` | Alturas dos painéis (px) | `780` / `560` / `520` / `760` |
| `REFERENCE_LINE_COLOR` / `COMPARISON_LINE_COLOR` | Cores das curvas | `"black"` / `"darkorange"` |
| `CLASS_ORDER` | Ordem das classes no diagnóstico | `["EN", "N", "LN"]` |
| `DECIMALS` | Casas decimais nas tabelas | `3` |

---

## ▶️ Como executar

**1. Calcular os índices e validar:**
```bash
python compute_roni_oni_psl_validation.py
```
Gera `ONI_timeseries.csv`, `RONI_timeseries.csv` e, se `RUN_VALIDATION=True`, os gráficos
comparativos com as métricas frente às séries oficiais.

**2. Comparar bases/diretórios de forma interativa:**
```bash
python compare_oni_roni_directories_plotly.py
```
Gera `ONI_RONI_directory_comparison.html`.

**3. Analisar skill, matriz de confusão e persistência:**
```bash
python generate_enso_cross_index_skill_persistence_html.py
```
Gera `ENSO_cross_index_skill_persistence.html`.

> Basta abrir os arquivos `.html` em qualquer navegador para explorar os gráficos interativos.

---

## 📤 Saídas

| Arquivo | Conteúdo |
|---|---|
| `ONI_timeseries.csv` | Colunas: `nino34_sst`, `nino34_anomaly`, `ONI`, `climatology_start`, `climatology_end` |
| `RONI_timeseries.csv` | Colunas: `nino34_anomaly`, `tropical_mean_anomaly`, `relative_nino34_raw`, `relative_nino34_scaled`, `RONI` |
| `ONI_RONI_directory_comparison.html` | Séries interativas e resumo consolidado de métricas entre bases |
| `ENSO_cross_index_skill_persistence.html` | Divergências, matriz de confusão, métricas de skill e persistência |

---

## 🧠 Metodologia (resumo)

- **ONI**: anomalia da Niño 3.4 medida contra uma **climatologia móvel de 30 anos** (atualizada a
  cada 5 anos, centrada nos blocos quinquenais), seguida de média móvel de 3 meses. Para o ano *y*,
  a janela é `[b(y)−15, b(y)+14]`, com `b(y) = 5·⌊(y−1)/5⌋ + 1`.
- **RONI**: anomalia da Niño 3.4 **menos** a anomalia média tropical (20°S–20°N), no mesmo instante,
  com ajuste de variância (`α = σ(N3.4)/σ(relativo)`) e média móvel de 3 meses. Remove boa parte da
  tendência global de aquecimento e reduz a dependência do período-base.
- **Classificação**: El Niño (índice ≥ +0,5 °C) / La Niña (≤ −0,5 °C) por ≥ 5 trimestres consecutivos;
  caso contrário, Neutro. Opcionalmente, estratificação por intensidade (fraco, moderado, forte,
  muito forte).
- **Validação**: RMSE, MAE, viés, correlação de Pearson (r), R², slope, intercept e desvio-padrão
  do erro contra as séries oficiais NOAA/PSL.

---

## 📚 Referências principais

- van Oldenborgh, G. J. et al. (2021). *Defining El Niño indices in a warming climate.* **Environ. Res. Lett.**, 16, 024003.
- L'Heureux, M. L. et al. (2024). *A Relative Sea Surface Temperature Index for Classifying ENSO Events in a Changing Climate.* **J. Climate**, 37(4), 1197–1211.
- Huang, B. et al. (2017). *Extended Reconstructed Sea Surface Temperature version 5 (ERSSTv5).* **J. Climate**, 30, 8179–8205.
- Trenberth, K. E. (1997). *The definition of El Niño.* **Bull. Amer. Meteor. Soc.**, 78(12), 2771–2777.
- Wilks, D. S. (2011). *Statistical Methods in the Atmospheric Sciences*, 3rd ed. Academic Press.
- Wilkinson, M. D. et al. (2016). *The FAIR Guiding Principles for scientific data management and stewardship.* **Scientific Data**, 3, 160018.
- NOAA/CPC — Páginas oficiais do ONI e do RONI (base ERSST, 1991–2020).

---

## 🔖 Citação

Se utilizar esta ferramenta, cite:

> Carbonel, A. (2026). *Uma ferramenta de código aberto para o cálculo e a avaliação dos índices
> ONI e RONI de monitoramento do El Niño–Oscilação Sul.*

---

## 📝 Licença

Sugestão: **MIT License** (adicione um arquivo `LICENSE` ao repositório).
Ajuste conforme a política de compartilhamento da sua instituição.

---

## 👤 Autora

**Alessandra Carbonel**
