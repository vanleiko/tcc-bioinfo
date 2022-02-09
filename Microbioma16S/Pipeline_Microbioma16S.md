## Microbioma 16S Workflow | Resumo

#### Bioinformática:
1. Baixar banco de dados.
2. Criar arquivos "metadata" e "manifest".
3. Importar sequências.
4. Análise da qualidade e inferência de amplicon.
5. Identificação taxonômica.
6. Resultado.

# Preparo do ambiente

- Instalar e ativar Qiime2:
```
wget https://data.qiime2.org/distro/core/qiime2-2021.8-py38-linux-conda.yml
conda env create -n qiime2-2021.8 --file qiime2-2021.8-py38-linux-conda.yml
rm qiime2-2021.8-py38-linux-conda.yml
```

# Download do banco de dados

- banco de dados SILVA:
```
wget https://data.qiime2.org/2021.8/common/silva-138-99-515-806-nb-classifier.qza
```

# Criar arquivos essenciais:

- Criar arquivo metadata:
```
echo -e "sample-id\tpaciente" > metadata-file.tsv
echo -e "#q2:types\tcategorical" >> metadata-file.tsv 
echo -e "MICROBIOMA\tjoao" >> metadata-file.tsv  
```

- Criar arquivo manifest:
```
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest-file.tsv
echo -e "MICROBIOMA\t$PWD/patient_joao_MICROBIOMA16S_S69_R1_001.fastq.gz\t$PWD/patient_joao_MICROBIOMA16S_S69_R2_001.fastq.gz" >> manifest-file.tsv
```

# Ativar ambiente
```
conda activate qiime2-2021.8
```

# Iniciar as análises
- Importar sequências:
```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest-file.tsv \
--output-path import.qza \
--input-format PairedEndFastqManifestPhred33V2 

qiime demux summarize --i-data import.qza --o-visualization import.qzv
```

- Controle de qualidade e inferência de Amplicon (sugestão: testar diferentes parâmetros após visualizar o “import.qzv”): 
```
qiime dada2 denoise-paired --i-demultiplexed-seqs import.qza \
 --p-trunc-len-f 250 \
 --p-trunc-len-r 250 \
 --p-trim-left-f 17 \
 --p-trim-left-r 21 \
 --o-representative-sequences rep-seqs.qza \
 --p-n-threads 20 \
 --o-table table.qza \
 --o-denoising-stats stats.qza
 
qiime metadata tabulate \
--m-input-file stats.qza \
--o-visualization stats.qzv
```

- Classificação toxomômica:
```
qiime feature-classifier classify-sklearn \
 --i-classifier silva-138-99-515-806-nb-classifier.qza\
 --i-reads rep-seqs.qza \
 --o-classification taxonomy.qza
 
qiime taxa barplot \
 --i-table table.qza \
 --i-taxonomy taxonomy.qza \
 --m-metadata-file metadata-file.tsv \
 --o-visualization taxa-bar-plots.qzv
```

# Desativar ambiente
```
conda deactivate
```

# Resultados:

### <b>Nível 1:</b>

Bacteria: 64,9%

Unassigned: 35,07%

Eucaryota: 0,03

![level1](https://user-images.githubusercontent.com/69684722/153259382-ea8659a4-57b4-4a1a-aec6-907c8a814cdf.png)

### <b>Nível 2: - Maiores frequências:</b>

Proteobacteria: 46,2%

Unassigned: 35,07%

Firmicutes: 18,05%

![level2](https://user-images.githubusercontent.com/69684722/153257398-87ee3c51-2d86-4059-be1f-b9ce139c0bc8.png)

### <b>Nível 6 - Maiores frequências:</b>

<i>Pseudomonas:<i/> 45,6%

Unassigned: 35,07%

<i>Enterococcus:</i> 13,87%

<i>Staphylococcus carnosus:</i> 3,67%

![level6](https://user-images.githubusercontent.com/69684722/153258298-a5a0cab8-c7b8-4fcd-8413-c6827ce90d8d.png)







 


