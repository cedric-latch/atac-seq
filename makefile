samples:
	basename data/raw/*_R1.fastq* | sed 's/_R1.fastq.*//g' > data/samples.txt
qc:
	mkdir -p data/qc && cat data/samples.txt |
	  parallel -j 4 fastqc data/raw/{}_R* -o data/qc 
trimgalore:
	mkdir -p data/trimgalore && cat data/samples.txt |
	  parallel -j 4 trim_galore --cores 1 --dont_gzip --output_dir data/trimgalore data/raw/{}_R*
index:
	mkdir -p data/bowtie2 && bowtie2-build data/GCF_000001405.26_GRCh38_genomic.fna hg38
align:
	ls data/samples.txt |
	  parallel bowtie2 \
		--very-sensitive \
		--threads 4 \
		-x data/bowtie2/hg38 \
		-U data/trimgalore/{}_R* \
		| samtools view -h -b - > data/bowtie2/{}.bam
