# Датасеты, бенчмарки и стратегия валидации

## Принцип валидации

SANG2-AMR должен быть валидирован по трёхуровневой схеме:

1. **Технический**: сигналы SANG2 корректно вычисляются на нуклеотидных последовательностях
2. **Биологический**: сигналы SANG2 имеют интерпретируемый биологический смысл
3. **Практический**: SANG2-AMR обнаруживает гены резистентности, которые пропускают reference-based инструменты

---

## Exploration-phase датасеты (первые 2 месяца)

### E1. Модельный организм — E. coli K-12 (MG1655)

**Зачем:** Самый изученный геном, ground truth для всех сигналов.
**Что проверяем:**
- MI half-life = 3 для кодирующих регионов (кодонная структура)
- ΔH показывает границы генов
- Autocorrelation lag-3 peak в ORF
- Regime classifier корректно разделяет кодирующие / интергенные

**Источник:**
- NCBI RefSeq: GCF_000005845.2
- GenBank: U00096.3
- Размер: 4.6 Mb, 4,288 генов с полной аннотацией

### E2. Плазмида — pCAV1392-131 (сигнальный контраст)

**Зачем:** Конъюгативная плазмида K. pneumoniae с транспозазами для контрастного анализа.
**Что проверяем:**
- Hurst / 3-mer entropy отличаются от хромосомы хозяина
- IS-элементы видны как signal gradient peaks
- Conjugal transfer region = coherent signal block

**Источник:**
- NCBI: CP011577.1 (pCAV1392-131, K. pneumoniae CAV1392)
- Размер: 130,719 bp, 15 транспозаз, 28 tra-генов

**NOTA BENE:** Ранее ошибочно указана как pKPC_CAV1321. Не содержит blaKPC.
Для E5 (AMR-детекция) нужна отдельная плазмида с blaKPC.

**Результаты exploration/002:** 6/7 сигналов различаются (p<0.01), transposase recall=100%.

### E3. Минимальный метагеном — HMP Mock Community

**Зачем:** Искусственное сообщество с известным составом (20 видов).
**Что проверяем:**
- Clustering разделяет контиги по видам без референса
- MI half-life стабилен внутри вида, отличается между видами

**Источник:**
- SRA: SRR2726667 (HMP even mock community)
- Ground truth: 20 видов с известными геномами

---

## Benchmark-phase датасеты (месяцы 2–4)

### B1. CARD (Comprehensive Antibiotic Resistance Database)

**Роль:** Ground truth для известных генов резистентности.
**Использование:** НЕ для обучения, а для валидации: «нашёл ли SANG2-AMR регион, где лежит известный ген?»
**Источник:** https://card.mcmaster.ca/download
**Размер:** ~5000 AMR генов, ~90,000 вариантов

### B2. ResFinder database

**Роль:** Второй independent ground truth (другая курация).
**Источник:** https://bitbucket.org/genomicepidemiology/resfinder_db
**Размер:** ~3,100 генов

### B3. CAMI Challenge datasets (metagenome benchmark)

**Зачем:** Стандартный бенчмарк для метагеномных инструментов.
**Что проверяем:** Clustering accuracy vs Kraken2, MetaPhlAn, VAMB.
**Источник:** https://data.cami-challenge.org/
**Релевантные:** CAMI II Marine, Strain Madness

### B4. PATRIC/BV-BRC AMR phenotype data

**Зачем:** Десятки тысяч изолятов с MIC-данными (Mode C validation).
**Что проверяем:** Корреляция SANG2-сигналов с фенотипом резистентности.
**Источник:** https://www.bv-brc.org/
**Размер:** >300,000 геномов с AMR-метаданными

---

## Валидационная матрица

| Эксперимент                          | Датасет      | Метрика                     | Baseline инструмент     |
|--------------------------------------|--------------|-----------------------------|--------------------------|
| MI half-life = 3 для CDS             | E1           | MAE от теоретического       | —                        |
| ΔH → границы генов                   | E1           | Precision/Recall vs аннот.  | Prodigal                 |
| Hurst различает хромосому и плазмиду | E1 + E2      | AUROC                       | PlasFlow, MOB-suite      |
| Cross-ref находит AMR-гены           | E2 + CARD    | Sensitivity/Specificity     | RGI, AMRFinderPlus       |
| Clustering метагенома                | E3, B3       | ARI (Adjusted Rand Index)   | VAMB, MetaBAT2           |
| Phenotype prediction (Mode C)        | B4           | AUROC по MIC                | Kover, Seq2Geno2Pheno    |
| **Novel AMR detection**              | Специальный* | Кол-во верифицированных     | RGI (negative control)   |

*Для novel detection: берём геномы с фенотипической резистентностью, где RGI не находит объяснения → SANG2-AMR должен предложить кандидатов.

---

## Стратегия «Novel AMR» валидации

Самая важная и самая сложная метрика — обнаружение **неизвестных** механизмов.

### Подход 1: Leave-one-out на CARD
1. Убираем один ген резистентности из CARD
2. Запускаем SANG2-AMR на геноме, содержащем этот ген
3. Проверяем: нашёл ли SANG2-AMR этот регион как аномальный?
4. Repeat для всех генов → Recall кривая

### Подход 2: Retrospective discovery
1. Берём геномы из PATRIC с «unexplained resistance» (фенотип R, но RGI ничего не нашёл)
2. Запускаем SANG2-AMR
3. Проверяем кандидатов вручную (BLAST, domain search, структурное предсказание)
4. Это уже не бенчмарк, а реальное открытие

### Подход 3: Synthetic insertion
1. Берём ген резистентности из одного вида
2. Вставляем в геном другого вида (in silico)
3. Меняем codon usage на хозяйский (усложнение)
4. Проверяем: находит ли SANG2-AMR вставку?
5. Варьируем степень «маскировки» → строим ROC-кривую

---

## Инструменты для сравнения (baselines)

| Задача                    | Baseline                     | Тип           |
|---------------------------|------------------------------|---------------|
| AMR gene detection        | RGI (CARD), AMRFinderPlus    | Reference-based |
| Metagenome binning        | VAMB, MetaBAT2, SemiBin     | ML + coverage  |
| Gene prediction           | Prodigal, Bakta              | Ab initio      |
| Plasmid detection         | PlasFlow, MOB-suite, gplas   | ML / MOB-typing |
| GI detection              | IslandViewer, Alien_Hunter   | Composition bias |
| Taxonomic classification  | Kraken2, MetaPhlAn4          | Reference-based |
| Phenotype prediction      | Kover, Seq2Geno2Pheno        | ML on k-mers   |

**Позиционирование SANG2-AMR:** Мы не конкурируем с reference-based инструментами.
Мы находим то, что они **пропускают**. Идеальная метрика:
`novel_candidates = SANG2_hits − (CARD ∪ ResFinder)`.
