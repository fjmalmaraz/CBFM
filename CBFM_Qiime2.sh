
### PRIMER PASO (ETAPA 1): generar un excel con las columnas 'sample-id', 'forward-absolute-filepath'
### y 'reverse-absolute-filepath' donde indicamos los ficheros y su localización, con '$' delante. p.ej:

###	Sample-id: FM019	
###	forward-absolute-filepath: $PWD/CBFM_Raw/21351-FM019-A-P214053-2-2389-16S_S29_L001_R1_001.fastq.gz
###	reverse-absolute-filepath: $PWD/CBFM_Raw/21351-FM019-A-P214053-2-2389-16S_S29_L001_R2_001.fastq.gz

### Para fusionar ficheros tenemos esta opción (sin necesidad de activar qiime2). 
### "DCIN" es la carpeta donde se encuentran los archivos originales.

#	mkdir Fusion
#	zless DCIN/FICHERO1.gz DCIN/FICHERO2.fastq.gz | 
#	gzip > Fusion/FIB038-A.fastq.gz


### Una vez tenemos el Excel hay que guardarlo en formato texto delimitado por tabulaciones (.txt), 
# al cual hay que forzarle el cambio a la terminación .tsv


### Con este fichero ejecutamos lo siguiente (primero activamos qiime escribiendo en la terminal: conda activate qiime2):

	qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 	--input-path CBFM_Raw_Manifest.tsv \
 	--input-format PairedEndFastqManifestPhred33V2 \
 	--output-path CBFM_01_sequences.qza

### Para ver la calidad de los ficheros, antes de hacer cutadapt o dada2 o lo que vayamos a hacer, lo visualizamos, para ello 
### generamos un fichero .qzv que luego visualizamos en qiime2.

	qiime demux summarize \
 	--i-data ./CBFM_01_sequences.qza \
  	--o-visualization ./CBFM_01_sequences.qzv

### Una vez generado el fichero .qza (en este caso, 'CBFM_01_sequences.qza') procedemos a eliminar el PRIMER de las lecturas,
### para ello existen dos formas, cutadapt y dada2 (No es exactamente así);
### Dada2 asigna a las lecturas un ASV que es una única Amplicon Sequence Variant una vez eliminado los ruidos. 
### Cutadapt elimina adaptadores, primers y Poly-A. 

### CUTADAPT:

	qiime cutadapt trim-paired \
 	--i-demultiplexed-sequences CBFM_01_sequences.qza \
  	--p-front-f CCTACGGGNGGCWGCAG \
  	--p-front-r GACTACHVGGGTATCTAATCC \
 	--p-minimum-length 200 \
 	--p-discard-untrimmed True \
 	--p-cores 20 \
	--o-trimmed-sequences CBFM_01_sequences_trimmed.qza

## FIJARSE BIEN QUE AHORA ESTAMOS COGIENDO LOS FICHEROS GENERADOS POR CUTADAPT! 

	qiime dada2 denoise-paired \
	--i-demultiplexed-seqs CBFM_01_sequences_trimmed.qza \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-trunc-len-f 250 \
	--p-trunc-len-r 250 \
	--p-max-ee-f 2 \
	--p-max-ee-r 2 \        #Expected errors
	--p-n-threads 0 \
	--o-table CBFM_01_table_with_cutadapt.qza \
	--o-representative-sequences CBFM_01_rep_seq_with_cutadapt.qza \
	--o-denoising-stats CBFM_01_denoising_with_cutadapt.qza


	qiime feature-table summarize \
	--i-table CBFM_01_table_with_cutadapt.qza \
	--o-visualization CBFM_01_table_with_cutadapt.qzv \
	--m-sample-metadata-file CBFM_Raw_Manifest.tsv


	qiime feature-table tabulate-seqs \
	--i-data CBFM_01_rep_seq_with_cutadapt.qza \
	--o-visualization CBFM_01_rep_seq_with_cutadapt.qzv


	qiime metadata tabulate \
	--m-input-file CBFM_01_denoising_with_cutadapt.qza \
	--o-visualization CBFM_01_denoising_with_cutadapt.qzv


## UNA VEZ QUE HEMOS LLEGADO A LOS FICHEROS CON CUTADAPT Y DADA2 Y COMPROBADO CALIDAD DE LOS FASTQ EN QIIME2 VIEW, PASAMOS A GENERAR LA BASE DE DATOS DE SILVA (https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494)

	DESDE LA TERMINAL:
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.txt.gz

gunzip tax_slv_ssu_138.1.txt.gz

	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/taxmap_slv_ssu_ref_nr_138.1.txt.gz

gunzip taxmap_slv_ssu_ref_nr_138.1.txt.gz

	
	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/tax_slv_ssu_138.1.tre.gz

gunzip tax_slv_ssu_138.1.tre.gz


	wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz

gunzip SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz

## AQUI AHORA ES SÚPER IMPORTANTE DESCOMPRIMIR LOS FICHEROS QUE ESTÁN EN FORMATO .gz ESCRIBIENDO EN LA TERMINAL: gunzip NOMBREFICHERO.gz

## CONTINUAMOS AHORA GENERANDO LAS BASES DE SILVA (ETAPA 2):
	
	qiime tools import \
    --type 'FeatureData[SILVATaxonomy]' \
    --input-path tax_slv_ssu_138.1.txt \
    --output-path taxranks-silva-138.1-ssu-nr99.qza 

	qiime tools import \
    --type 'FeatureData[SILVATaxidMap]' \
    --input-path taxmap_slv_ssu_ref_nr_138.1.txt \
    --output-path taxmap-silva-138.1-ssu-nr99.qza 

	qiime tools import \
    --type 'Phylogeny[Rooted]' \
    --input-path tax_slv_ssu_138.1.tre \
    --output-path taxtree-silva-138.1-nr99.qza 

	qiime tools import \
    --type 'FeatureData[RNASequence]' \
    --input-path SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta \
    --output-path silva-138.1-ssu-nr99-rna-seqs.qza

	qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

	
	qiime rescript parse-silva-taxonomy \
    --i-taxonomy-tree taxtree-silva-138.1-nr99.qza \
    --i-taxonomy-map taxmap-silva-138.1-ssu-nr99.qza \
    --i-taxonomy-ranks taxranks-silva-138.1-ssu-nr99.qza \
    --o-taxonomy silva-138.1-ssu-nr99-tax.qza

## Seguimos en SILVA, ahora estamos con el proceso: “Culling” low-quality sequences with cull-seqs;

	qiime rescript cull-seqs \
    --i-sequences silva-138.1-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

## Seguimos en SILVA, ahora el proceso es: Filtering sequences by length and taxonomy	

	qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.1-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

## Seguimos en SILVA, ahora el proceso es: Dereplicating in uniq mode

	qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.1-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

# Seguimos en SILVA; Cambiado los primers conforme a los datos de nuestros primers. Creemos que que el primero es la 341F (337F cuando busco en Escholeria Echoli). Reverse primer encontrado en 789R.(GGATTAGATACCCTGGTAGTC) Secuencia complementaria y reversa 

	qiime feature-classifier extract-reads \
    --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer CCTACGGGNGGCWGCAG \
    --p-r-primer GACTACHVGGGTATCTAATCC \
    --p-n-jobs 25 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.1-ssu-nr99-seqs-515f-806r.qza

	qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-515f-806r.qza \
    --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza \
    --o-dereplicated-taxa  silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza			## Sí, este paso se hace dos comandos antes, pero
                                                                                        #si nos fijamos bien, estamos cogiendo los ficheros
                                                                                        #generados por el comando anterior justamente. Javi
                                                                                        #también lo procesa 2 veces, por eso yo también.
	qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.1-ssu-nr99-seqs-515f-806r-uniq.qza \
    --i-reference-taxonomy silva-138.1-ssu-nr99-tax-515f-806r-derep-uniq.qza \
    --o-classifier silva-138.1-ssu-nr99-515f-806r-classifier.qza


## Una vez generada la base de SILVA; continuamos ahora con la asignación taxonómica (ETAPA 3);

	qiime feature-classifier   classify-sklearn \
    --i-reads cutadaptwithdada2_24124/CBFM_01_rep_seq_with_cutadapt.qza \
    --i-classifier SILVA24125/silva-138.1-ssu-nr99-515f-806r-classifier.qza \
    --p-n-jobs 20 \
    --o-classification TAXONOMYTERESA.qza

### Para generar una tabla con la clasificación

	qiime metadata tabulate \
      --m-input-file TAXONOMYTERESA.qza \   ##
 	 --o-visualization TAXONOMYTERESA.qzv

	##Esto es para meterle el excel con las variables, pero de momento no se las voy a meter yo.
 	--m-input-file METADATA \ 


### Permite visualizar los diferentes grupos taxonómicos por participante en la muestra y usar qiime 2 para ver los gráficos
qiime taxa barplot \
  --i-table cutadaptwithdada2_24124/CBFM_01_table_with_cutadapt.qza \
  --i-taxonomy TAXONOMYTERESA.qza \
  --o-visualization CBFM-taxa-bar-plotsTERESA.qzv

	--m-input-file METADATA		##Este METADATA es el de los datos extras para cruzar variables

### Agrupar por muestras al nivel de género (5)

qiime taxa collapse \
  --i-table cutadaptwithdada2_24124/CBFM_01_table_with_cutadapt.qza \
  --i-taxonomy TAXONOMYTERESA.qza \
  --p-level 6 \
  --o-collapsed-table  CBFM-genusTERESA.qza

qiime metadata tabulate \
      --m-input-file CBFM-genusTERESA.qza \
      --o-visualization CBFM-genusTERESA.qzv

	## --m-input-file m-input-file METADATATERESA \

qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences cutadaptwithdada2_24124/CBFM_01_rep_seq_with_cutadapt.qza \
      --o-alignment CBFM-aligned-rep-seqsTERESA.qza \
      --o-masked-alignment CBFM-mask-aligned-red-seqsTERESA.qza \
      --o-tree unrooted-treeTERESA.qza \
      --o-rooted-tree rooted-treeTERESA.qza



