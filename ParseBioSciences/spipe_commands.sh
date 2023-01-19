# spipe.v1.0.0 is already installed as conda env
conda activate spipe.v1.0.0

PBS='/data1/ivanir/Ilaria2023/ParseBS'

PATH="/data1/ivanir/Ilaria2023/ParseBS/ParseBiosciences-Pipeline.1.0.3p:$PATH"

# Genome indexing
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.gtf.gz \
--output_dir $PBS/newvolume/genomes/hg38

# Pipeline running 
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/S1 

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_2.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_2.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/S2 

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-22604.UnspecifiedIndex.HNLNNDRX2.s_2.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-22604.UnspecifiedIndex.HNLNNDRX2.s_2.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/nuclei 
