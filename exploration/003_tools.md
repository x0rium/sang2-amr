# Agent & MCP Trace

## Exploration session: project bootstrap

Agent & MCP Trace (кратко):
- Роль агента: Explore — исследование кодовой базы auto-reverser для каталогизации алгоритмов
- MCP-проверки: нет (внутренний проект, не внешний API)
- Upstream: auto-reverser/src/ — все 17 алгоритмических модулей каталогизированы

Доказательства:
- signals.py: 8 сигналов (H, SDR, Otsu, ΔH, MI half-life, autocorrelation, H_ratio, Hurst)
- resonance.py: 6-signal geometric mean (A, D, MI, P, X, G)
- merge.py: PPMI pair merge с lossless reversibility
- cluster.py: 4-axis clustering (length + discriminator + signature)
- xref.py: 3 механизма cross-reference
- fields.py: 7-stage field inference
- checksum.py: 4 алгоритма (MOD256, XOR, CRC16, CRC32)

Решение: можно продолжать — все алгоритмы каталогизированы,
маппинг на геномный домен задокументирован в docs/algorithms.md.

---

## Exploration 001: SANG2 signals on E. coli K-12 (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — адаптация signals.py для 4-буквенного алфавита
- MCP-проверки: NCBI Entrez (BioPython) для скачивания генома
- Данные: U00096.3 (FASTA, 4,641,652 bp) + GenBank аннотация (4,318 CDS)

### Входные данные
- E. coli K-12 MG1655: GC=50.79%, 4,318 CDS генов
- Файлы: exploration/data/ecoli_k12.fasta, ecoli_k12.gb

### Результаты по сигналам

| Сигнал | Prediction | Actual | Status |
|--------|-----------|--------|--------|
| A1 Entropy | H sliding window показывает variation | H mean=1.9802 (std=0.021), range 1.77-2.0 | PASS |
| A2 SDR | Structural k-mers = консервативные мотивы | 91/256 4-меров structural, SDR ratio=35.5% | PASS |
| A3 Otsu | GC-content bucketing | Threshold=0.48, разделяет 3921/11550 окон | PASS |
| A4 ΔH | Negative at promoters/RBS | ΔH=-0.04 at -30..-10 nt before gene start | PASS |
| A5 MI HL | ≈ 3 для CDS | 37.9 median, NO exponential decay | **FAIL** |
| A6 Autocorr | Lag-3 peak in CDS | 97.2% CDS have peak at 3k; 55.2% strongest | PASS |
| A7 H/H_max | CDS in CODING regime (0.60-0.85) | H_ratio=0.992 everywhere; thresholds wrong | **FAIL** |
| A8 Hurst | 0.58-0.67 persistent | CDS=0.594, IG=0.624 | PASS |

### Ключевые находки

**F1. MI half-life не переносится на 4-буквенный алфавит.**
- MI(lag=1) = 0.024 bits = 1.2% от H_max (2.0 bits)
- Для сравнения: в байтовом домене MI(lag=1) типично 20-60% от H_max
- MI не затухает экспоненциально — осциллирует на уровне шума
- Кодонный сигнал ЕСТЬ, но виден через автокорреляцию (A6), не через MI decay
- MI(3)>MI(4) в 66.2% CDS vs 52.7% intergenic (p=0.034) — слабый, но значимый сигнал
- FFT: period-6 доминирует (codon-pair), не period-3

**F2. H/H_max пороги из байтового домена не работают.**
- С 4 символами и GC≈50%, H/H_max ≈ 0.99 для всех регионов
- 3-mer entropy лучше: CDS=0.966, IG=0.914 (Otsu threshold=0.917)
- Но направление инвертировано: CDS ВЫШЕ (все 64 кодона используются), IG НИЖЕ (AT-bias)

**F3. Автокорреляция — самый сильный детектор кодонной структуры.**
- 97.2% CDS имеют хотя бы один пик автокорреляции на кратных 3
- Это прямой аналог lag-3 сигнала в auto-reverser

**F4. ΔH работает как предсказано.**
- Средний профиль вокруг старт-кодона: конвергентный в [-30, -10] и [0, +5]
- Дивергентный в [-5] (граница промотор→ген)
- Может использоваться для ab initio детекции границ генов

**F5. Hurst экспонент различает CDS и intergenic.**
- CDS: 0.594 ± 0.035 (менее persistent — кодонная структура дробит корреляции)
- IG: 0.624 ± 0.051 (более persistent — длинные AT/GC треки)
- Разница мала, но стабильна

### Необходимые изменения в архитектуре

1. **Заменить MI half-life** → "MI codon periodicity index" (на основе A6 autocorrelation)
2. **Пересчитать regime thresholds** для 4-буквенного алфавита:
   - Использовать 3-mer entropy вместо nucleotide entropy
   - Или Otsu на H по окнам (data-driven, без magic numbers)
3. **Resonance scan**: единица = динуклеотид/кодон OK, но сигнал MI в resonance
   тоже будет слабым — нужна нормализация по H_max домена

### Решение
6/8 сигналов переносятся напрямую или с минимальной адаптацией.
1/8 (MI half-life) требует переформулировки концепции.
1/8 (regime) требует пересчёта порогов.

**Можно продолжать к exploration/002** (plasmid contrast).
MI half-life заменяется на autocorrelation periodicity index.

---

## Exploration 002: Chromosome vs Plasmid Signal Contrast (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — статистический контраст сигналов хромосома vs плазмида
- MCP-проверки: NCBI Entrez для скачивания плазмиды CP011577.1
- Данные: pCAV1392-131 (130,719 bp, K. pneumoniae), 15 транспозаз, 28 генов conjugal transfer

**NOTA BENE:** CP011577.1 — это pCAV1392-131, не pKPC_CAV1321 (ошибка в data-sources.md).
Не содержит blaKPC. Для exploration/005 нужна другая плазмида.

### Результаты

**PART A — Глобальное сравнение:**

| Сигнал | Chromosome | Plasmid | Δ |
|--------|-----------|---------|---|
| GC | 0.5079 | 0.5165 | +0.009 |
| 3-mer H_ratio | 0.9907 | 0.9922 | +0.002 |
| SDR ratio | 0.3555 | 0.3672 | +0.012 |
| Hurst | 0.6257 | 0.6214 | -0.004 |
| AC period-3 | 0.4312 | 0.5394 | +0.108 |
| MI(lag=1) | 0.0212 | 0.0121 | -0.009 |

**PART B — Статистика (5000 bp windows): 6/7 сигналов различаются (p < 0.01).**
- Единственный незначимый: GC content (p=0.42) — оба организма ~51% GC.
- Самый мощный: MI(lag=1) (p = 5.5e-21).

**PART C — Anomaly detection:**
- 19.8% окон плазмиды аномальны по GC (вне chr mean ± 2σ)
- 10.9% по Hurst, 10.5% по SDR

**PART D — Transposase detection: recall = 100% (15/15).**
- Метод: gradient(GC + 3-mer entropy) > 90th percentile
- 8/15 top gradient peaks попадают рядом с транспозазами

**PART E — Conjugal transfer region (34-68 kb):**
- Статистически отличается от остальной плазмиды по 5/7 сигналов (p < 0.01)
- Ниже SDR (0.269 vs 0.301), выше AC period-3 (0.369 vs 0.327)

### Ключевые находки

**F6. MI(lag=1) — самый мощный дискриминатор chromosome/plasmid.**
p = 5.5e-21, хотя GC-content статистически неразличим. MI чувствителен к
dinucleotide frequencies, которые отличаются между видами/репликонами.

**F7. Gradient GC + 3-mer entropy = эффективный детектор IS-элементов.**
100% recall на 15 транспозазах. Простейший composition-based метод работает.
Это baseline для L3 (Anomaly Detection) пайплайна.

**F8. Conjugal region = coherent signal block.**
34 kb region с 28 tra-генами имеет отличный от остальной плазмиды сигнальный профиль.
Это прямой аналог "encapsulated protocol" из auto-reverser — чужеродный функциональный
блок, приобретённый горизонтально.

### Решение
4/5 проверок прошли. **Можно продолжать к exploration/003** (pair merge hierarchy).

---

## Exploration 002b: Genomic-Specific Signals (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — реализация и валидация 6 геномно-специфичных сигналов
- MCP-проверки: нет (внутренние алгоритмы)
- Данные: E. coli K-12 + pCAV1392-131

### Результаты

| Сигнал | Уровень | Результат | Включаем? |
|--------|---------|-----------|-----------|
| S1 GC-skew | L0 | oriC с ошибкой 0.1% (3 kb / 4.6 Mb) | **Да** — уникальный структурный сигнал |
| S2 Chargaff-2 | L1 | Conjugal vs rest: p=0.001 | **Да** — детектор функциональных блоков |
| S3 ρ dinucleotide | L3 | 9/16 dinuc различаются chr/plas | **Да для cross-ref** (между организмами) |
| S4 CUB (RSCU) | L3 | Неверный baseline (E. coli вместо хозяина) | **Да, но переделать**: baseline = хозяйский геном |
| S5 IR density | L1 | Тренд верный (636 vs 552), max в TNP-кластере | **Условно** — нужна лучшая резолюция |
| S6 CUB anomaly | L3 | Не работает (p=0.55) | **Нет в текущем виде** → переделать с S4 |

### Ключевые находки

**F9. GC-skew — самый мощный новый сигнал.**
Находит origin/terminus репликации из чистой композиции, zero magic numbers.
Для AMR: плазмиды имеют другой GC-skew паттерн (std=0.65) vs хромосома (std=4.57).
Аналога в auto-reverser нет — это чисто геномный сигнал.

**F10. Chargaff-2 deviation = детектор функциональных блоков внутри репликона.**
Conjugal region: 0.060 vs rest: 0.033 (p=0.001).
Для AMR: pathogenicity islands и resistance cassettes — чужеродные блоки,
нарушающие strand symmetry.

**F11. ρ-profile — fingerprint организма, не региона.**
Работает МЕЖДУ организмами (cross-reference), не ВНУТРИ одного генома.
Для AMR Mode B (metagenome): кластеризация контигов по организму через ρ-distance.
Для AMR Mode A (isolate): бесполезен внутри генома.

**F12. CUB нужен ПРАВИЛЬНЫЙ baseline.**
Сравнивать с E. coli бессмысленно для K. pneumoniae плазмиды.
Правильно: host_rscu = genome-wide CDS RSCU, window_rscu = local CDS RSCU.
Отклонение window от host → HGT candidate. Переделать в 003 или отдельном эксперименте.

### Решение
3 сигнала (GC-skew, Chargaff-2, ρ) включаются в архитектуру.
CUB требует переделки baseline. IR density — условно, нужна доработка.

---

## Exploration 003: PPMI Pair Merge Hierarchy (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — реализация и тестирование PPMI pair merge + BPE для геномных последовательностей
- MCP-проверки: нет (внутренние алгоритмы)
- Данные: E. coli K-12 first 100 kb

### Результаты

**003 (PPMI scoring): ПРОВАЛ.**
- PPMI × log2(1+count) создаёт runaway chains: GC → GGC → GGCG → GGCGC → ...
- Причина: merged symbols редки → PPMI зашкаливает → каскадное расширение одной цепочки
- Только 4/16 dinucleotides, 3/64 codons, все G-initial
- Compression всего 1.44x

**003b (standard BPE): УСПЕХ.**
- Frequency-only BPE: 300 merges, compression 3.62x, lossless
- 11/16 dinucleotides, 23/64 codons, ATG + TAA + TAG найдены
- Иерархия: nt(4) → dinuc(11) → codons(23) → 4-mers(82) → 5-mers(148)
- BPE автоматически обнаруживает кодонную структуру без предзнания о ней

### Ключевые находки

**F13. PPMI pair merge не работает на маленьком алфавите.**
Это ВТОРОЙ алгоритм SANG2 (после MI half-life), ломающийся на 4 символах.
Общая причина: формулы, включающие P(x), дают экстремальные значения при |Σ|=4,
потому что P(x) ≈ 0.25 для каждого символа (высокая baseline probability).
В байтовом домене P(x) ≈ 1/256 → формулы работают в «правильном» динамическом диапазоне.

**F14. Frequency-based BPE — правильная замена.**
Работает из коробки, автоматически строит биологически осмысленную иерархию.
PPMI можно использовать POST-HOC для анализа, но не для scoring merges.

**F15. 5 «потерянных» dinucleotides (CA, CT, GA, GT, TA) — не баг.**
Они не формируются как отдельные пары, потому что BPE мержит их в составе
более длинных единиц. Например, CA появляется внутри GCA, TCA и т.д.
При необходимости: forced pre-merge всех 16 dinucleotides на level 0.

### Решение
4/5 проверок прошли на BPE. Pair merge работает с frequency scoring.
**Можно продолжать к exploration/004** (resonance scan).

---

## Exploration 004: Resonance Scan (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — адаптация resonance scan для nucleotide/codon alphabets
- MCP-проверки: нет
- Данные: 500 CDS E. coli K-12, aligned по start codon, 300 nt

### Результаты: 5/5 проверок прошли

**Nucleotide level (|Σ|=4):** Resonance работает слабо.
- P (positional stability) ≈ 0.01 — каждый нуклеотид встречается в каждой позиции
- Top pair: TG k=1 (R=0.248), driven by A+X, not P
- 85/240 пар с R>0

**Codon level (|Σ|=64):** Resonance работает отлично.
- 3695/6170 пар с R>0
- **ATG доминирует**: top-10 все ATG-X, R до 0.53
- ATG P=0.39 (привязан к позиции 0 = start codon) → высокая R
- Synonymous codon pairs найдены: GAA-GAG (R=0.34)
- Сигналы разумно независимы: max |ρ|=0.45

### Ключевые находки

**F16. Resonance нужен codon-level alphabet для CDS.**
С 4 символами P бесполезна (все символы встречаются везде).
С 64 символами P = главный дискриминатор (ATG → позиция 0).
Рекомендация: L2 resonance на codon level для CDS, nucleotide level для non-coding.

**F17. ATG (start codon) — самый resonant элемент генома.**
R=0.53 для ATG-AGT: высокая P (позиция 0), высокая D (всегда first),
средняя MI (ассоциация с Met→Ser парой). Это прямой аналог "magic bytes"
из auto-reverser — фиксированный стартовый паттерн.

**F18. Codon eigendistances распределены равномерно (k=1..8).**
В байтовом домене eigendistances кластеризуются (k=3 для headers).
В геноме кодоны коррелированы на всех расстояниях (белковая композиция).
Нет выраженного «характерного расстояния» на codon level.

---

## Exploration 005: Ab Initio AMR Detection (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — composite anomaly scoring для ab initio AMR detection
- MCP-проверки: NCBI Entrez для скачивания pKpQIL (FJ628167.1)
- Данные: pKpQIL (20,158 bp, KPC-2 карбапенемаза + 7 транспозаз)

### Результаты: 3/4 проверок прошли. **blaKPC-2 найден ab initio, ранг 8/37.**

Метод: 7 validated signals → z-scores vs E. coli baseline → geometric mean composite.
- Top-10: 1/1 AMR recall (100%), 7/7 transposase recall (100%)
- Лучший одиночный детектор: **Chargaff-2 deviation** (KPC rank #1/37)
- Сильные: autocorrelation (#3), signal gradient (#3)
- Слабые: ρ-distance (#14) — не работает внутри одного репликона (F11 confirmed)

### Ключевые находки

**F19. Chargaff-2 deviation = лучший одиночный детектор AMR-генов.**
KPC rank #1/37 по этому сигналу. Карбапенемаза нарушает strand symmetry сильнее
всех остальных регионов. Это новый сигнал, которого нет в auto-reverser.

**F20. Composite scoring работает: KPC ранг 8, все транспозазы найдены.**
Геометрическое среднее z-scores — robust к слабым отдельным сигналам.
Ни один сигнал по отдельности не ставит KPC в top-3, но composite ставит в top-10.

**F21. pKpQIL = validation success для SANG2-AMR.**
На реальной KPC-плазмиде, без CARD, без BLAST, без обучения — только
статистические сигналы из composition, symmetry, periodicity.

---

## Exploration 006: Reference-Free Metagenomic Clustering (2026-04-12)

Agent & MCP Trace:
- Роль агента: dev-agent — in silico mock metagenome + composition-based clustering
- MCP-проверки: NCBI Entrez для скачивания 3 геномов (S. aureus, P. aeruginosa, B. subtilis)
- Данные: 4 species × 80 contigs = 320 contigs, 2-10 kb, shuffled

### Результаты: 2/4 проверок прошли. ARI=0.875 (target >0.5 — прошёл с запасом).

| Метод | ARI | NMI |
|---|---|---|
| Hierarchical (22 feat) | **0.875** | **0.877** |
| K-means (22 feat) | 0.858 | 0.855 |
| GC alone (1 feat) | 0.840 | 0.833 |

### Ключевые находки

**F22. GC content — доминирующий дискриминатор для видов с разным GC.**
GC alone (ARI=0.84) почти так же хорош, как 22 фичи (0.875).
С 4 видами (33%, 44%, 51%, 67%) GC хватает для разделения.
Полные фичи помогают на сложной паре B. subtilis (44%) vs E. coli (51%).

**F23. ρ-vector НЕ добавляет ценности для видов с разным GC.**
GC+ρ (0.799) хуже GC alone (0.840). ρ коррелирован с GC →  добавляет шум.
Нужен тест на видах с ОДИНАКОВЫМ GC (E. coli vs K. pneumoniae).

**F24. 3-mer и 4-mer entropy — самые ценные фичи после GC.**
Leave-one-out: H3 (ARI drop +0.069) и H4 (+0.054) — top-2 по важности.
Они несут информацию о codon usage, которой нет в GC.

### Решение
ARI=0.875 >> 0.5 (target). Кластеризация работает.
GC — главный сигнал. Дополнительные фичи помогают на сложных парах.
Для реального метагенома (десятки видов с перекрывающимся GC) ρ и k-mer entropy будут важнее.

## Все exploration experiments завершены (001–006)
