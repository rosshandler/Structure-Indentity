# spipe.v1.0.0 is already installed as conda env
conda activate spipe.v1.0.0

PBS='/data1/ivanir/Ilaria2023/ParseBS'

PATH="/data1/ivanir/Ilaria2023/ParseBS/ParseBiosciences-Pipeline.1.0.3p:$PATH"

cd $PBS/newvolume/genomes/
cat Homo_sapiens.GRCh38.108.gtf EmGFP.gtf > Homo_sapiens.GRCh38.108.EmGFP.gtf
cat Homo_sapiens.GRCh38.dna.primary_assembly.fa EmGFP.fa > Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa

# Genome indexing
split-pipe \
--mode mkref \
--genome_name hg38 \
--fasta $PBS/newvolume/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.EmGFP.fa \
--genes $PBS/newvolume/genomes/Homo_sapiens.GRCh38.108.EmGFP.gtf \
--output_dir $PBS/newvolume/genomes/hg38

# Merge both lanes of single cell experiment
cd $PBS/newvolume/expdata/
cat SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_1.r_1.fq.gz SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_2.r_2.fq.gz > SLX-22601.r_1.fq.gz
cat SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_1.r_1.fq.gz SLX-22601.UnspecifiedIndex.HGTNNDMXY.s_2.r_2.fq.gz > SLX-22601.r_2.fq.gz

nohup /data2/ivanir/Feline2023/ParseBS/newvolume/expdata/demultiplexer.rhel/demuxFQ \
    -c -d -i -e -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22601.lostreads.r_1.fq.gz \
    -s SLX-22601.demultiplexsummary.r1.txt \
    SLX-22601.r_1.index.txt \
    SLX-22601.r_1.fq.gz &
    
nohup /data2/ivanir/Feline2023/ParseBS/newvolume/expdata/demultiplexer.rhel/demuxFQ \
    -c -d -i -e -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22601.lostreads.r_2.fq.gz \
    -s SLX-22601.demultiplexsummary.r2.txt \
    SLX-22601.r_2.index.txt \
    SLX-22601.r_2.fq.gz &

rm SLX-22601.r_1.fq.gz SLX-22601.r_2.fq.gz

# Pipeline running
#single cell
split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.AGTCAA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.AGTTCC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.ATGTCA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.CTTGTA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.GATCAG.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.GATCAG.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

split-pipe --mode all --kit WT --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFast/SLX-22601.DNAA007.HGMLNDMXY.TAGCTT.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sCell

#single nuclei
nohup /home/ivanir/tools/demultiplexer.rhel/demultiplexer.rhel/demuxFQ \
    -c -d -i -e -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22604.lostreads.r_1.fq.gz \
    -s SLX-22604.demultiplexsummary.r1.txt \
    SLX-22604.r_1.index.txt \
    SLX-22604.UnspecifiedIndex.HNLNNDRX2.s_2.r_1.fq.gz &
    
nohup /home/ivanir/tools/demultiplexer.rhel/demultiplexer.rhel/demuxFQ \
    -c -d -i -e -t 1 -r 0.01 \
    -o correctFastq \
    -b SLX-22604.lostreads.r_2.fq.gz \
    -s SLX-22604.demultiplexsummary.r2.txt \
    SLX-22604.r_2.index.txt \
    SLX-22604.UnspecifiedIndex.HNLNNDRX2.s_2.r_2.fq.gz &

split-pipe --mode all --kit WT_mini --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22604.HGMLNDMXY.ACTTGA.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22604.HGMLNDMXY.ACTTGA.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sNuclei/ACTTGA

split-pipe --mode all --kit WT_mini --chemistry v2 --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/correctFastq/SLX-22604.HGMLNDMXY.CAGATC.s_1.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/correctFastq/SLX-22604.HGMLNDMXY.CAGATC.s_1.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/sNuclei/CAGATC



