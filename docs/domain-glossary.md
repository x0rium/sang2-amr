# Глоссарий AMR / Геномика → SANG2

Для разработчика, приходящего из мира протоколов и reverse engineering.

---

## Молекулярная биология (базовые понятия)

**Нуклеотид (nt)** — элементарный символ ДНК: A, C, G, T. Аналог: один байт.

**Кодон** — тройка нуклеотидов, кодирующая одну аминокислоту. 64 кодона → 20 аминокислот + 3 стоп-кодона. Аналог: 3-байтовый opcode с lookup-таблицей.

**Codon bias** — неравномерное использование синонимичных кодонов. Каждый организм предпочитает свои варианты. Аналог: разные реализации протокола используют разные encoding для одного и того же поля.

**ORF (Open Reading Frame)** — участок ДНК от стартового кодона (ATG) до стоп-кодона. Кандидат в ген. Аналог: сообщение протокола от header до footer.

**Промотор** — регуляторная последовательность перед геном (-10 box: TATAAT, -35 box: TTGACA). Аналог: magic bytes в начале пакета.

**RBS (Ribosome Binding Site)** — Shine-Dalgarno последовательность (AGGAGG), сигнал начала трансляции. Аналог: sync byte / preamble.

**Терминатор** — последовательность, завершающая транскрипцию (часто палиндром → шпилька). Аналог: end-of-message delimiter.

---

## Мобильные генетические элементы

**Плазмида** — кольцевая ДНК, реплицирующаяся независимо от хромосомы. Часто несёт гены резистентности. Аналог: отдельный протокол, работающий параллельно с основным.

**Транспозон** — «прыгающий ген», перемещающийся внутри и между геномами. Окружён inverted repeats (IR). Аналог: переиспользуемый блок кода, копирующийся между проектами.

**IS-элемент (Insertion Sequence)** — простейший транспозон: transposase gene + flanking IR. Аналог: самоперемещающийся payload с magic bytes на концах.

**Интегрон** — система захвата генов. Содержит integrase + attI site + gene cassettes. Кассеты вставляются/удаляются. Аналог: extensible протокол с TLV-подобными плагинами.

**Gene cassette** — мобильная единица: один ген + attC site (59-base element). Аналог: TLV-запись (type-length-value).

**Pathogenicity Island (PAI)** — крупный блок (10-200 kb) чужеродной ДНК с генами вирулентности/резистентности. Характерные признаки: другой GC%, tRNA на фланге, прямые повторы. Аналог: encapsulated protocol — чужой протокол внутри основного потока.

---

## Антимикробная резистентность

**MIC (Minimum Inhibitory Concentration)** — минимальная концентрация антибиотика, подавляющая рост. Чем выше MIC → тем резистентнее. Аналог: порог отказа — при каком давлении протокол перестаёт работать.

**β-лактамаза** — фермент, разрушающий β-лактамные антибиотики (пенициллины, цефалоспорины, карбапенемы). Наиболее распространённый механизм резистентности. Семейства: TEM, SHV, CTX-M, KPC, NDM, OXA.

**ESBL (Extended-Spectrum β-Lactamase)** — β-лактамаза с расширенным спектром. Устойчивость к цефалоспоринам 3-го поколения.

**Карбапенемаза** — β-лактамаза, разрушающая карбапенемы (антибиотики «последней линии»). KPC, NDM, VIM, OXA-48. Критическая угроза по ВОЗ.

**MRSA (Methicillin-Resistant S. aureus)** — стафилококк с геном mecA (изменённый PBP2a). Парадигма AMR.

**VRE (Vancomycin-Resistant Enterococcus)** — энтерококк с кластером van генов. Изменяет мишень антибиотика.

**mcr (mobilized colistin resistance)** — ген резистентности к колистину (последний резервный антибиотик). Обнаружен в 2015, на плазмидах. Пример «нового механизма, которого не было в базах».

---

## Базы данных и инструменты

**CARD (Comprehensive Antibiotic Resistance Database)** — курируемая база AMR-генов + онтология (ARO). ~5000 генов. Инструмент: RGI (Resistance Gene Identifier).

**ResFinder** — датская база AMR-генов (CGE, DTU). Фокус на приобретённой резистентности.

**AMRFinderPlus** — инструмент NCBI. Использует HMM-профили + BLAST.

**NCBI SRA (Sequence Read Archive)** — крупнейший архив сырых секвенирований.

**BV-BRC (Bacterial and Viral Bioinformatics Resource Center)** — бывший PATRIC. Геномы + фенотипы AMR.

**CAMI (Critical Assessment of Metagenome Interpretation)** — бенчмарк для метагеномных инструментов. Community challenge.

---

## Секвенирование

**Illumina** — short-read секвенирование. Риды 150–300 bp, высокая точность (>99.9%), высокая пропускная способность. Стандарт для клинических лабораторий.

**Oxford Nanopore (ONT)** — long-read секвенирование. Риды до 1 Mb (!), точность 95–99%, реалтайм. Идеален для замыкания плазмид и разрешения повторов.

**PacBio HiFi** — long-read с высокой точностью (>99.9%). Риды 10–25 kb. Золотой стандарт для полных геномов.

**Coverage / глубина покрытия** — сколько раз каждая позиция генома прочитана. Типично 30–100x для isolate, 1–10x для metagenome. Плазмиды обычно имеют **более высокий** coverage (больше копий на клетку).

---

## Маппинг терминологии для быстрого перехода

| SANG2 / auto-reverser    | AMR / Геномика                     |
|--------------------------|-------------------------------------|
| Alphabet: byte (256)     | Alphabet: nucleotide (4) or aa (20) |
| Packet                   | Read / Contig                       |
| Session                  | Sample / Isolate                    |
| Message type             | Gene family / Operon                |
| Field                    | Domain / Motif                      |
| Magic bytes              | Promoter / RBS / Start codon        |
| Delimiter                | Stop codon / Terminator             |
| Checksum                 | Conserved flanking regions          |
| Cross-reference          | Horizontal Gene Transfer            |
| Encapsulated protocol    | Pathogenicity Island                |
| TLV record               | Gene cassette (in integron)         |
| Protocol fuzzing         | In silico mutagenesis               |
| Unknown protocol         | Novel resistance mechanism          |
