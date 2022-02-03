## Metagenômica Workflow | Resumo

#### Bioinformática:
1. Análise da qualidade dos dados.
2. Remoção de contaminação do hospedeiro humano.
3. Identificação do patógeno.
4. Relatório final.

# Preparo do Ambiente

- Download miniconda:

```wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh```

- Instalar miniconda:
 
```bash Miniconda3-py38_4.10.3-Linux-x86_64.sh```

- Instalar dependências do ambiente:

```conda install --channel defaults conda python=3.6 --yes```

```conda update --channel defaults --all --yes```

- Instalar os programas:
 
```conda install --channel bioconda fastqc cutadapt kraken2 krona bwa samtools spades --yes```

- Criar diretórios para organização dos arquivos:

```mkdir -p fastq kraken2 kraken2-db fastqc cutadapt human-ref bwa```

# Download dos dados

- Banco de dados Kraken2:

```wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz```

- Genoma humano de referência hg38 UCSC:

```wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz```

# Reorganizar arquivos em seus respectivos diretórios

- Mover os arquivos Fastq:

```mv *fastq.gz fastq/ ```

- Descompactar o banco de dados Kraken2 e movê-lo para seu respectivo diretório:

```tar -xf k2_standard_8gb_20210517.tar.gz --directory kraken2-db/```

- Mover o genoma humano de referência para seu respectivo diretório:

```mv Homo_sapiens_UCSC_hg38.tar.gz human-ref/```

# Controle de qualidade e limpeza das sequências

- Gerar relatórios de qualidade do sequenciamento com fastqc:

```fastqc fastq/METAGENOMICA_R1_001.fastq.gz fastq/METAGENOMICA_R2_001.fastq.gz -o fastqc/```

- Filtragem e trimagem das sequências com cutadapt:

```cutadapt -u 9 -U 9 -u -1 -U -1 -m 50 \```

```-o cutadapt/METAGENOMICA_cleaned_R1.fastq.gz -p cutadapt/METAGENOMICA_cleaned_R2.fastq.gz \```

```fastq/METAGENOMICA_R1_001.fastq.gz fastq/METAGENOMICA_R2_001.fastq.gz > cutadapt/summary_cutadapt_METAGENOMICA.txt```

# Remover contaminantes do hospedeiro

- Proceder com o mapeamento contra o genoma humano de referência (converter SAM para BAM direto - não precisa salvar o SAM):

```bwa mem human-ref/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa \```

```cutadapt/METAGENOMICA_cleaned_R1.fastq.gz cutadapt/METAGENOMICA_cleaned_R2.fastq.gz \```

```| samtools view -b > bwa/METAGENOMICA_mapped_host.bam ```

- Gerar BAM dos reads não mapeados:

```samtools view -u -f 12 -b bwa/METAGENOMICA_mapped_host.bam \```

```| samtools sort -n > bwa/METAGENOMICA_unmapped_host.bam```

- Gerar relatório com os counts dos reads mapeados em humano:

```samtools flagstat bwa/METAGENOMICA_mapped_host.bam > bwa/METAGENOMICA_mapped_host_flagstat.txt```

- Gerar FASTq dos reads não mapeados:

```samtools fastq bwa/METAGENOMICA_unmapped_host.bam \```

```-1 bwa/METAGENOMICA_unmapped_host_R1.fastq -2 bwa/METAGENOMICA_unmapped_host_R2.fastq```

# Identificação taxonômica

- Identificação taxonômica dos reads:

```kraken2 -db kraken2-db/ \```

```--report kraken2/METAGENOMICA_kraken2_report.txt \```

```--output kraken2/METAGENOMICA_kraken2_NT.out \```

```--minimum-base-quality 20 \```
  
```--paired bwa/METAGENOMICA_unmapped_host_R1.fastq bwa/METAGENOMICA_unmapped_host_R2.fastq```

- Gerar relatório:

```ktImportTaxonomy kraken2/METAGENOMICA_kraken2_NT.out \```

```-q 2 -t 3 -o kraken2/METAGENOMICA_classification.html```

# Montagem dos reads em contigs (etapa opcional de validação)

- Criar diretório para os outputs do SPAdes:

```mkdir assembly_metagenomica```

- Montagem dos contigs:

```nohup spades.py --meta \```
 
 ```-1 bwa/METAGENOMICA_unmapped_host_R1.fastq \```
  
 ```-2 bwa/METAGENOMICA_unmapped_host_R2.fastq \```
  
 ```-o assembly_metagenomica & ```
 
- Identificação taxonômica:

```kraken2 -db kraken2-db/ \```

```--report kraken2/METAGENOMICA_kraken2_report_contigs.txt \```

```--output kraken2/METAGENOMICA_kraken2_NT_contigs.out \```

```--minimum-base-quality 20 \```

```assembly_metagenomica/contigs.fasta```

- Gerar relatório:

```ktImportTaxonomy kraken2/METAGENOMICA_kraken2_NT_contigs.out \```

```-q 2 -t 3 -o kraken2/METAGENOMICA_classification_contigs.html```

# Inspecionar resultados

Nesta última etapa, cabe ao analista ou médico patologista/infectologista inspecionar os resultados gerados pelos relatórios de diversidade para identificar possíveis patógenos na amostra.
É importante ressaltar que o resultado da identificação taxonômica gerado pelo Kraken2 não deve ser levado em consideração unicamente. O processo de laudamente deve envolver uma etapa posterior de validação do achado, que poderia remover falso-positivos gerados pelo Kraken2 ou mesmo encontrar patógenos que não foram apontado pela identificação do Kraken2.

## Classificação taxonômica - reads:

![snapshot](https://user-images.githubusercontent.com/69684722/152216945-3b7cdec1-e828-489d-a79d-784b7e459915.svg)

## Classificação taxonômica - contigs:

![snapshot (1)](https://user-images.githubusercontent.com/69684722/152217001-677fc9ac-9198-466c-a616-9fd2f94ba65d.svg)

