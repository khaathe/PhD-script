### Step 1: Building the STAR index.*

STAR --runMode genomeGenerate --genomeDir <star_index_path> --genomeFastaFiles <reference> --sjdbOverhang 100 --sjdbGTFfile <gencode.v22.annotation.gtf> --runThreadN 8

### Step 2: Alignment 1st Pass.

STAR --genomeDir <star_index_path> --readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,... --runThreadN <runThreadN> --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --readFilesCommand <bzcat|cat|zcat> --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None

### Step 3: Intermediate Index Generation.

STAR --runMode genomeGenerate --genomeDir <output_path> --genomeFastaFiles <reference> --sjdbOverhang 100 --runThreadN <runThreadN> --sjdbFileChrStartEnd <SJ.out.tab from previous step>

### Step 4: Alignment 2nd Pass.

STAR --genomeDir <output_path from previous step> --readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,... --runThreadN <runThreadN> --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --readFilesCommand <bzcat|cat|zcat> --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline <formatted RG line provided by wrapper>

htseq-count -m intersection-nonempty -i gene_id -r pos -s no  ./GNS/G144.bam ../../GDC_Ref/gencode.v22.annotation.gtf