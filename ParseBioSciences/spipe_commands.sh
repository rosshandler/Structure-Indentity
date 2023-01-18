# spipe is already installed as conda env
conda activate spipe


PBS='/data1/ivanir/Chp2022/ParseBS'

# Genome indexing --> see Chp repository

# Pipeline running 
split-pipe --mode all --kit WT_mini --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-21921.DNAA007.H52HVDRX2.s_2.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-21921.DNAA007.H52HVDRX2.s_2.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/DNAA008 

split-pipe --mode all --kit WT_mini --genome_dir $PBS/newvolume/genomes/hg38/ \
--fq1 $PBS/newvolume/expdata/SLX-21921.DNAA008.H52HVDRX2.s_2.r_1.fq.gz \
--fq2 $PBS/newvolume/expdata/SLX-21921.DNAA008.H52HVDRX2.s_2.r_2.fq.gz \
--output_dir $PBS/newvolume/analysis/DNAA008 
