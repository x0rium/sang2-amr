# Архитектура SANG2-AMR

## Обзор пайплайна

```
NGS-риды (FASTQ)
    │
    ▼
┌──────────────────────────────────────────────────────┐
│  L0: Signal Extraction                               │
│  SANG2-universal (на скользящем окне по контигам):   │
│    Энтропия (nt + 3-mer), SDR, ΔH, автокорреляция,  │
│    Hurst, MI(lag=1)                                  │
│  Genomic-specific:                                    │
│    GC-skew (кумулятивный), Chargaff-2 deviation,     │
│    ρ dinucleotide odds ratio                         │
│  [MI half-life НЕ переносится — см. gotcha G15]      │
│  [H/H_max regime: использовать 3-mer H — см. G16]   │
└─────────────┬────────────────────────────────────────┘
              │
              ▼
┌──────────────────────────────────────────────────────┐
│  L1: Structure Detection                             │
│  SANG2-universal:                                     │
│    Pair merge → иерархия единиц                      │
│    Relative peak segmentation → атомы                │
│    Regime classification (3-mer entropy) → зонирование│
│  Genomic-specific:                                    │
│    Chargaff-2 блоки → функциональные регионы         │
│    Inverted repeat density → IS-element boundaries   │
└─────────────┬────────────────────────────────────────┘
              │
              ▼
┌──────────────────────────────────────────────────────┐
│  L2: Resonance Scan (universal, 6 signals)           │
│  A D MI P X G pair scan по нуклеотидам/кодонам       │
│  → мотивы, eigendistances, домены                    │
│  → кандидаты в функциональные единицы                │
│  [Resonance остаётся alphabet-agnostic]              │
└─────────────┬────────────────────────────────────────┘
              │
              ▼
┌──────────────────────────────────────────────────────┐
│  L3: Anomaly Detection (AMR-specific)                │
│  SANG2-universal:                                     │
│    Cross-reference → мобильные элементы              │
│    Hurst anomaly → чужеродные вставки                │
│    Regime shift → pathogenicity islands               │
│  Genomic-specific:                                    │
│    ρ-distance from host → HGT candidates             │
│    CUB anomaly (RSCU vs host) → foreign genes        │
│    Signal gradient (GC + 3-mer H) → IS-elements      │
└─────────────┬────────────────────────────────────────┘
              │
              ▼
┌──────────────────────────────────────────────────────┐
│  L4: Clustering & Classification                     │
│  Length bucketing + discriminators                    │
│  ρ-profile clustering (Mode B: metagenome binning)   │
│  → группировка по функции/организму                  │
│  → сопоставление с known AMR (опц.)                  │
└─────────────┬────────────────────────────────────────┘
              │
              ▼
┌──────────────────────────────────────────────────────┐
│  L5: Report                                          │
│  Кандидаты в новые AMR-механизмы                     │
│  Мобильные элементы + контекст                       │
│  Визуализация сигналов                               │
└──────────────────────────────────────────────────────┘
```

**Принцип разделения:** Resonance (L2) остаётся **universal** — работает на любом алфавите
с 6 pair-сигналами. Геномно-специфичные знания (strand symmetry, codon usage,
replication geometry) живут в L0/L1/L3 и не загрязняют универсальное ядро.

## Входные данные

### Первичный вход
- **FASTQ** — сырые риды (Illumina short-read, Nanopore/PacBio long-read)
- **FASTA** — собранные контиги/скаффолды (после SPAdes/MEGAHIT/Flye)
- **BAM/SAM** — выровненные риды (опционально, для позиционного анализа)

### Метаданные (опционально)
- MIC (Minimum Inhibitory Concentration) — фенотип резистентности
- Источник образца (клинический, окружающая среда, животноводство)
- Антибиотик, к которому тестируется резистентность

## Предобработка (вне ядра SANG2)

Не реализуем сами — используем стандартные инструменты:

1. **QC**: fastp / Trimmomatic → обрезка адаптеров, фильтрация по качеству
2. **Сборка**: SPAdes (short-read) / Flye (long-read) / MEGAHIT (metagenome)
3. **Предсказание ORF**: Prodigal / Pyrodigal → координаты генов-кандидатов
4. **Baseline AMR**: RGI/CARD или AMRFinderPlus → известные гены (для сравнения)

SANG2-AMR получает на вход **контиги + ORF-аннотацию** и работает с ними.

## Режимы работы

### Mode A: Single Isolate
Один бактериальный геном. Цель: найти нестандартные гены резистентности,
мобильные элементы, pathogenicity islands.

### Mode B: Metagenome
Метагеномная проба (кишечник, почва, сточные воды). Цель: профиль резистома
без полной таксономической классификации. Clustering на L4 группирует
фрагменты по организму-источнику.

### Mode C: Cohort (панель изолятов)
N изолятов с MIC-данными. Цель: найти геномные детерминанты, коррелирующие
с фенотипом резистентности (GWAS-like, но на уровне структурных сигналов).

## Интеграция с auto-reverser

Ядро SANG2 (signals, resonance, merge, cluster) переиспользуется напрямую.
Адаптационный слой:

```
auto-reverser/src/          amr/src/
├── signals.py         →    genome_signals.py    (алфавит A/C/G/T, окно = контиг)
├── resonance.py       →    genome_resonance.py  (пары нуклеотидов/кодонов)
├── merge.py           →    genome_merge.py      (иерархия: nt→codon→motif→domain→gene)
├── cluster.py         →    genome_cluster.py    (кластеризация ридов/контигов)
├── xref.py            →    hgt_detect.py        (горизонтальный перенос = cross-reference)
├── fsm.py             →    не используется напрямую
├── fields.py          →    region_classify.py   (классификация геномных регионов)
├── checksum.py        →    conserved_flanks.py  (консервативные фланки = "чексуммы" генов)
└── stream.py          →    contig_parser.py     (парсинг FASTA/FASTQ → внутренний IR)
```

## IR (Internal Representation) для геномики

```python
@dataclass
class Contig:
    id: str
    sequence: bytes          # A=0, C=1, G=2, T=3 (4-буквенный алфавит)
    quality: bytes | None    # Phred scores, если из FASTQ
    source: str              # sample ID

@dataclass
class GenomicRegion:
    contig_id: str
    start: int
    end: int
    strand: Literal['+', '-']
    kind: str                # 'orf', 'intergenic', 'repeat', 'mobile_element', 'unknown'
    signals: SignalVector    # все 8 SANG2-сигналов для этого региона

@dataclass
class SignalVector:
    entropy: float
    sdr_ratio: float
    delta_h: float
    mi_halflife: float
    autocorr_peaks: list[int]
    h_ratio: float
    regime: str              # SPARSE | RAW_SIGNAL | SYMBOLIC | HIGH_ENTROPY
    hurst: float

@dataclass
class AMRCandidate:
    region: GenomicRegion
    evidence: list[str]      # какие сигналы сработали и почему
    resonance_score: float   # R из geometric mean
    mobility_score: float    # вероятность мобильного элемента
    novelty: str             # 'known' | 'variant' | 'novel'
    nearest_known: str | None  # ближайший известный ген из CARD (если есть)
    confidence: float        # 0-1
```

## Принципы архитектуры

1. **Данные → пороги.** Ни одного magic number. Всё через SDR/Otsu.
2. **Add-on, не замена.** Стандартные инструменты на preprocessing, SANG2 на анализ.
3. **Один проход.** Все L0-сигналы вычисляются за один проход по контигу (vectorized).
4. **Lossless hierarchy.** Pair merge сохраняет split path — можно вернуться к нуклеотидам.
5. **Resonance = relational.** Минимальная единица для resonance — пара (dinucleotide / codon pair). Никогда single nucleotide.
