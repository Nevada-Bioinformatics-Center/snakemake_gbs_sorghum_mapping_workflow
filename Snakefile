import glob

configfile: "config.yaml"

inputdirectory=config["input_data"]
ref=config["ref"]
print(inputdirectory)
SAMPLES1, SAMPLES =glob_wildcards(inputdirectory+"/Sample_{sample1}/{sample}_R1_clipped_passed-re-filter.fastq.bz2", followlinks=True)

print(SAMPLES)

##### target rules #####
rule all:
    input: 
       "calls/all_merged.vcf.gz",
       "calls/all_merged.q20dp8.vcf.gz",
       directory("plots/origfile"),
       directory("plots/file1"),
       directory("plots/file2"),
       #"calls/all_merged.q20dp8saf0rpr1.hmp.txt",
       #"calls/all_merged.q20dp8.hmp.txt"
       "calls/all_merged.q20dp8.ped",
       "calls/all_merged.q20dp8.imiss",
       "calls/all_merged.q20dp8.012",
       "calls/all_merged.q20dp8saf0rpr1.ped",
       "calls/all_merged.q20dp8saf0rpr1.imiss",
       "calls/all_merged.q20dp8mis50.imiss",
       "calls/all_merged.q20dp8mis25.imiss",
       "calls/all_merged.noBadSamples.maf05dp8mis50.imiss",
       "calls/all_merged.q20dp8mis75.imiss",
       "calls/all_merged.q20dp8saf0rpr1.012",
       "calls/all_merged.noBadSamples.maf05dp8mis50.012",
       "calls/all_merged.noBadSamples.maf05dp8mis50.ped",
       "all_samples.txt",
       "qc/multiqc_report_all.html"


rule bwa_index:
    input:
        ref
    output:
        idx=multiext("{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/genome.log"
    params:
        prefix="genome",
        algorithm="bwtsw"
    resources: time_min=520, mem_mb=20000, cpus=1
    wrapper:
        "v1.7.0/bio/bwa/index"

rule uncompress_fastq_r1:
    input:
        inputdirectory+"/Sample_{sample}/{sample}_R1_clipped_passed-re-filter.fastq.bz2"
    output:
        temp("uncompressed/{sample}_R1_clipped_passed-re-filter.fastq")
    log:
        "logs/uncompress/{sample}_R1.log"
    threads: 1
    shell:
        "bzcat {input} > {output} 2> {log}"

rule uncompress_fastq_r2:
    input:
        inputdirectory+"/Sample_{sample}/{sample}_R2_clipped_passed-re-filter.fastq.bz2"
    output:
        temp("uncompressed/{sample}_R2_clipped_passed-re-filter.fastq")
    log:
        "logs/uncompress/{sample}_R2.log"
    threads: 1
    shell:
        "bzcat {input} > {output} 2> {log}"

rule fastqc_posttrim_r1:
    input:
        "uncompressed/{sample}_R1_clipped_passed-re-filter.fastq"
    output:
        html="qc/fastqc_posttrim/{sample}_r1.html",
        zip="qc/fastqc_posttrim/{sample}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{sample}_r1.log"
    resources: time_min=520, mem_mb=20000, cpus=1
    threads: 1
    wrapper:
        "0.73.0/bio/fastqc"

rule fastqc_posttrim_r2:
    input:
        "uncompressed/{sample}_R2_clipped_passed-re-filter.fastq"
    output:
        html="qc/fastqc_posttrim/{sample}_r2.html",
        zip="qc/fastqc_posttrim/{sample}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{sample}_r2.log"
    threads: 1
    resources: time_min=520, mem_mb=20000, cpus=1
    wrapper:
        "0.73.0/bio/fastqc"


rule bwa_mem:
    input:
        #reads=[inputdirectory+"/Sample_{sample}/{sample}_R1_clipped_passed-re-filter.fastq.bz2", inputdirectory+"/Sample_{sample}/{sample}_R2_clipped_passed-re-filter.fastq.bz2"]
        reads=[ancient("uncompressed/{sample}_R1_clipped_passed-re-filter.fastq"), ancient("uncompressed/{sample}_R2_clipped_passed-re-filter.fastq")],
        idx=multiext("genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "mapped/{sample}.sorted.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index="genome",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 16
    resources: time_min=1320, mem_mb=20000, cpus=16
    wrapper:
        "v1.7.0/bio/bwa/mem"


rule samtools_index:
    input:
        ancient("mapped/{sample}.sorted.bam")
    output:
        "mapped/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.73.0/bio/samtools/index"


rule sambamba_merge:
    input:
       ancient(expand("mapped/{sample}.sorted.bam", sample=SAMPLES))
    output:
        "mapped/all_merged.bam"
    log:
        "logs/sambamba-merge/all_merge.log"
    params:
        "" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        16     # This value - 1 will be sent to -@
    resources: time_min=1320, mem_mb=20000, cpus=16
    wrapper:
        "0.74.0/bio/sambamba/merge"

rule samtools_index_merged:
    input:
        ancient("mapped/all_merged.bam")
    output:
        "mapped/all_merged.bam.bai"
    params:
        "" # optional params string
    threads: 1
    resources: time_min=520, mem_mb=20000, cpus=1
    wrapper:
        "0.73.0/bio/samtools/index"

rule freebayes:
    input:
        ref=ref,
        # you can have a list of samples here
        samples=ancient("mapped/all_merged.bam"),
        # the matching BAI indexes have to present for freebayes
        indexes=ancient("mapped/all_merged.bam.bai")
        # optional BED file specifying chromosomal regions on which freebayes 
        # should run, e.g. all regions that show coverage
        #regions="/path/to/region-file.bed"
    output:
        "calls/all_merged.vcf"  # either .vcf or .bcf
    log:
        "logs/freebayes/all_merged.log"
    params:
        #extra="--min-base-quality 10 --min-supporting-allele-qsum 10 --read-mismatch-limit 3 --min-coverage 5 --min-alternate-count 4 --exclude-unobserved-genotypes --genotype-qualities --ploidy 2 --no-mnps --no-complex --no-indels --mismatch-base-quality-threshold 10",         # optional parameters
        extra="--min-base-quality 10 --min-supporting-allele-qsum 10 --read-mismatch-limit 3 --min-coverage 5 --min-alternate-count 4 --exclude-unobserved-genotypes --genotype-qualities --mismatch-base-quality-threshold 10 --ploidy 2 -X -i -u",         # optional parameters
        chunksize=100000, # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # flag to use bcftools norm to normalize indels
    resources: time_min=4320, mem_mb=40000, cpus=32
    threads: 32
    wrapper:
        "0.73.0/bio/freebayes"


rule bgzip_vcf:
    input:
        "calls/all_merged.vcf"
    output:
        "calls/all_merged.vcf.gz"
    threads: 1
    resources: time_min=220, mem_mb=4000, cpus=1
    log:
        "logs/bgzip/bgzip.log"
    conda:
        "envs/htslib.yaml"
    shell:
        "bgzip {input} {output} 2> {log}"

rule filter_vcf:
    input:
        "calls/all_merged.vcf.gz"
    output:
        "calls/all_merged.q20dp8.vcf.gz"
    log:
        "logs/filter/vcftools_q20dp8.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    params:
        #Filters for minquality of 20, 
        #min depth of 8 (from LGC), minimum allele frequency of 0.05
        extra="--minQ 20 --minDP 8 --maf 0.05 --recode-INFO-all"
    wrapper:
        "0.74.0/bio/vcftools/filter"


rule filter_vcf2:
    input:
        "calls/all_merged.vcf.gz"
    output:
        "calls/all_merged.q20dp8mis50.vcf.gz"
    log:
        "logs/filter/vcftools_q20dp8miss50.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    params:
        #Filters for minquality of 20, 
        #min depth of 8 (from LGC), minimum allele frequency of 0.05
        extra="--minQ 20 --minDP 8 --maf 0.05 --max-missing 0.5 --recode-INFO-all"
    wrapper:
        "0.74.0/bio/vcftools/filter"

rule filter_vcf3:
    input:
        "calls/all_merged.vcf.gz"
    output:
        "calls/all_merged.q20dp8mis25.vcf.gz"
    log:
        "logs/filter/vcftools_q20dp8miss50.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    params:
        #Filters for minquality of 20, 
        #min depth of 8 (from LGC), minimum allele frequency of 0.05
        extra="--minQ 20 --minDP 8 --maf 0.05 --max-missing 0.25 --recode-INFO-all"
    wrapper:
        "0.74.0/bio/vcftools/filter"


rule filter_vcf4:
    input:
        "calls/all_merged.vcf.gz"
    output:
        "calls/all_merged.q20dp8mis75.vcf.gz"
    log:
        "logs/filter/vcftools_q20dp8miss75.log"
    resources: time_min=220, mem_mb=10000, cpus=1
    threads: 1
    params:
        #Filters for minquality of 20, 
        #min depth of 8 (from LGC), minimum allele frequency of 0.05
        extra="--minQ 20 --minDP 8 --maf 0.05 --max-missing 0.75 --recode-INFO-all"
    wrapper:
        "0.74.0/bio/vcftools/filter"

rule filter_vcf_strands:
    input:
        "calls/all_merged.q20dp8.vcf.gz"
    output:
        "calls/all_merged.q20dp8saf0rpr1.vcf.gz"
    log:
        "logs/filter/vcflib_strands.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcflib.yaml"
    shell:
        # SAF > 0 SAR > 0 filters for  reads on both strands
        # RPR > 1 RPL > 1 filters for at least two reads "balanced" on each side of the site
        "zcat {input} | vcffilter -f \"SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1\"  | bgzip > {output} 2> {log}"

rule bcftstats_origfile:
    input:
        vcf="calls/all_merged.vcf.gz",
        fa=ref
    output:
        "stats/all_merged.vcf.gz.stats"
    log:
        "logs/filter/bcfstats_origfile.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats -F {input.fa} -s - {input.vcf} > {output} 2> {log}"


rule bcftstats_file1:
    input:
        vcf="calls/all_merged.q20dp8.vcf.gz",
        fa=ref
    output:
        "stats/all_merged.q20dp8.vcf.gz.stats"
    log:
        "logs/filter/bcfstats_file1.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats -F {input.fa} -s - {input.vcf} > {output} 2> {log}"

rule bcftstats_file2:
    input:
        vcf="calls/all_merged.q20dp8saf0rpr1.vcf.gz",
        fa=ref
    output:
        "stats/all_merged.q20dp8saf0rpr1.vcf.gz.stats"
    log:
        "logs/filter/bcfstats_file2.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats -F {input.fa} -s - {input.vcf} > {output} 2> {log}"

rule bcfplot_origfile:
    input:
        "stats/all_merged.vcf.gz.stats"
    output:
        directory("plots/origfile")
    log:
        "logs/filter/bcfplot_origfile.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "plot-vcfstats -p {output} {input} 2> {log}"

rule bcfplot_file1:
    input:
        "stats/all_merged.q20dp8.vcf.gz.stats"
    output:
        directory("plots/file1")
    log:
        "logs/filter/bcfplot_file1.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "plot-vcfstats -p {output} {input} 2> {log}"

rule bcfplot_file2:
    input:
        "stats/all_merged.q20dp8saf0rpr1.vcf.gz.stats"
    output:
        directory("plots/file2")
    log:
        "logs/filter/bcfplot_file2.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "plot-vcfstats -p {output} {input} 2> {log}"

rule vcf_to_plink1:
    input:
        "calls/all_merged.q20dp8.vcf.gz"
    output:
        "calls/all_merged.q20dp8.ped",
        "calls/all_merged.q20dp8.map"
    log:
        "logs/vcftools/convert_file1.log"
    params:
        prefix="calls/all_merged.q20dp8"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --plink --out {params.prefix} 2> {log}"

rule vcf_to_plink2:
    input:
        "calls/all_merged.q20dp8saf0rpr1.vcf.gz"
    output:
        "calls/all_merged.q20dp8saf0rpr1.ped",
        "calls/all_merged.q20dp8saf0rpr1.map"
    log:
        "logs/vcftools/convert_file2.log"
    params:
        prefix="calls/all_merged.q20dp8saf0rpr1"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --plink --out {params.prefix} 2> {log}"

rule vcf_to_numeric1:
    input:
        "calls/all_merged.q20dp8.vcf.gz"
    output:
        "calls/all_merged.q20dp8.012",
    log:
        "logs/vcftools/convert_filenum1.log"
    params:
        prefix="calls/all_merged.q20dp8"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --012 --out {params.prefix} 2> {log}"

rule vcf_to_numeric2:
    input:
        "calls/all_merged.q20dp8saf0rpr1.vcf.gz"
    output:
        "calls/all_merged.q20dp8saf0rpr1.012",
    log:
        "logs/vcftools/convert_file2.log"
    params:
        prefix="calls/all_merged.q20dp8saf0rpr1"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --012 --out {params.prefix} 2> {log}"

rule vcf_missingindtest:
    input:
        "calls/all_merged.q20dp8mis50.vcf.gz"
    output:
        "calls/all_merged.q20dp8mis50.imiss",
    log:
        "logs/vcftools/missingtest.log"
    params:
        prefix="calls/all_merged.q20dp8mis50"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.prefix} 2> {log}"

rule vcf_miss50:
    input:
        "calls/all_merged.q20dp8mis50.imiss",
    output:
        "calls/lowDP.indv"
    log:
        "logs/lowdp.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/mawk.yaml"
    shell:
        "mawk '$5 > 0.7' {input} | cut -f1 > {output} 2> {log}"

rule filter_bybadseqsamp:
    input:
        vcf="calls/all_merged.vcf.gz",
        highmiss="calls/lowDP.indv"
    output:
        "calls/all_merged.noBadSamples.vcf.gz"
    log:
        "logs/filter/vcftools_noBadSamples.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input.vcf} --remove {input.highmiss} --recode --recode-INFO-all --stdout | gzip > {output}"


rule filter_bybadseqsampmis50:
    input:
        "calls/all_merged.noBadSamples.vcf.gz"
    output:
        "calls/all_merged.noBadSamples.maf05dp8mis50.vcf.gz"
    log:
        "logs/filter/vcftools_noBadSamplesMaxMiss50.log"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    params:
        #Filters for minquality of 20, 
        #min depth of 8 (from LGC), minimum allele frequency of 0.05
        extra="--minQ 20 --minDP 8 --maf 0.05 --max-missing 0.5 --recode-INFO-all"
    wrapper:
        "0.74.0/bio/vcftools/filter"

rule vcf_missingindnobad:
    input:
        "calls/all_merged.noBadSamples.maf05dp8mis50.vcf.gz"
    output:
        "calls/all_merged.noBadSamples.maf05dp8mis50.imiss"
    log:
        "logs/vcftools/missingtest_nobadsamp.log"
    params:
        prefix="calls/all_merged.noBadSamples.maf05dp8mis50"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.prefix} 2> {log}"

rule vcf_to_numericnobad:
    input:
        "calls/all_merged.noBadSamples.maf05dp8mis50.vcf.gz"
    output:
        "calls/all_merged.noBadSamples.maf05dp8mis50.012"
    log:
        "logs/vcftools/convert_filenumeric_nobadsamples.log"
    params:
        prefix="calls/all_merged.noBadSamples.maf05dp8mis50"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --012 --out {params.prefix} 2> {log}"

rule vcf_to_plinknobad:
    input:
        "calls/all_merged.noBadSamples.maf05dp8mis50.vcf.gz"
    output:
        "calls/all_merged.noBadSamples.maf05dp8mis50.ped"
    log:
        "logs/vcftools/convert_fileplinknobad.log"
    params:
        prefix="calls/all_merged.noBadSamples.maf05dp8mis50"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --plink --out {params.prefix} 2> {log}"


rule vcf_missingind3:
    input:
        "calls/all_merged.q20dp8mis25.vcf.gz"
    output:
        "calls/all_merged.q20dp8mis25.imiss",
    log:
        "logs/vcftools/missingtest.log"
    params:
        prefix="calls/all_merged.q20dp8mis25"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.prefix} 2> {log}"

rule vcf_missingind4:
    input:
        "calls/all_merged.q20dp8mis75.vcf.gz"
    output:
        "calls/all_merged.q20dp8mis75.imiss",
    log:
        "logs/vcftools/missingtest75.log"
    params:
        prefix="calls/all_merged.q20dp8mis75"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.prefix} 2> {log}"

rule vcf_missingind1:
    input:
        "calls/all_merged.q20dp8.vcf.gz"
    output:
        "calls/all_merged.q20dp8.imiss",
    log:
        "logs/vcftools/missing1.log"
    params:
        prefix="calls/all_merged.q20dp8"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.prefix} 2> {log}"

rule vcf_missingind2:
    input:
        "calls/all_merged.q20dp8saf0rpr1.vcf.gz"
    output:
        "calls/all_merged.q20dp8saf0rpr1.imiss",
    log:
        "logs/vcftools/missing2.log"
    params:
        prefix="calls/all_merged.q20dp8saf0rpr1"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.prefix} 2> {log}"

rule multiqc:
    input:
        expand("qc/fastqc_posttrim/{sample}_r1_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc_posttrim/{sample}_r2_fastqc.zip", sample=SAMPLES),
        expand("mapped/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "qc/multiqc_report_all.html"
    log:
        "logs/multiqc_all.log"
    resources: time_min=520, mem_mb=40000, cpus=1
    wrapper:
        "0.73.0/bio/multiqc"

rule make_samplelist:
    input:
        vcf="calls/all_merged.noBadSamples.maf05dp8mis50.vcf.gz"
    output:
        "all_samples.txt"
    resources: time_min=220, mem_mb=8000, cpus=1
    threads: 1
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools query -l {input.vcf} | sort > {output}"

#rule make_samplelist:
#    params: samples = SAMPLES
#    output:
#        "all_samples.txt",
#    resources: time_min=20, mem_mb=800, cpus=1
#    threads: 1
#    script:
#       "scripts/print_samples.py"
#


#rule sort1:
#    input:
#        "calls/all_merged.q20dp8.vcf.gz"
#    output:
#        "calls/all_merged.q20dp8.sorted.vcf.gz"
#    log:
#        "logs/tassel/sort_file1.log"
#    resources: time_min=220, mem_mb=24000, cpus=1
#    params:
#        javamem="22"
#    threads: 1
#    conda:
#        #"envs/tassel.yaml"
#        "envs/bcftools.yaml"
#    shell:
#        #"run_pipeline.pl -Xmx{params.javamem}G -SortGenotypeFilePlugin -inputFile {input} -outputFile {output} -fileType VCF 2> {log}"
#        "bcftools sort -m {params.javamem}G -Oz {input} -o {output} 2> {log}"
#
#rule sort2:
#    input:
#        "calls/all_merged.q20dp8saf0rpr1.vcf.gz"
#    output:
#        "calls/all_merged.q20dp8saf0rpr1.sorted.vcf.gz"
#    log:
#        "logs/tassel/sort_file2.log"
#    resources: time_min=220, mem_mb=24000, cpus=1
#    params:
#        javamem="22"
#    threads: 1
#    conda:
#        "envs/bcftools.yaml"
#    shell:
#        #"run_pipeline.pl -Xmx{params.javamem}G -SortGenotypeFilePlugin -inputFile {input} -outputFile {output} -fileType VCF 2> {log}"
#        "bcftools sort -m {params.javamem}G -Oz {input} -o {output} 2> {log}"
#
#rule tassel_convert1:
#    input:
#        "calls/all_merged.q20dp8.sorted.vcf.gz"
#    output:
#        "calls/all_merged.q20dp8.hmp.txt"
#    log:
#        "logs/tassel/convert_file1.log"
#    resources: time_min=220, mem_mb=24000, cpus=1
#    params:
#        prefix="calls/all_merged.q20dp8",
#        javamem="22"
#    threads: 1
#    conda:
#        "envs/tassel.yaml"
#    shell:
#        "run_pipeline.pl -Xmx{params.javamem}G -FileLoadPlugin -sortPositions true -vcf {input} -export {output} -exportType Hapmap 2> {log}"
#
#rule tassel_convert2:
#    input:
#        "calls/all_merged.q20dp8saf0rpr1.sorted.vcf.gz"
#    output:
#        "calls/all_merged.q20dp8saf0rpr1.hmp.txt"
#    log:
#        "logs/tassel/convert_file2.log"
#    resources: time_min=220, mem_mb=24000, cpus=1
#    params:
#        prefix="calls/all_merged.q20dp8saf0rpr1",
#        javamem="22"
#    threads: 1
#    conda:
#        "envs/tassel.yaml"
#    shell:
#        "run_pipeline.pl -Xmx{params.javamem}G -FileLoadPlugin -sortPositions true -vcf {input} -export {output} -exportType Hapmap 2> {log}"

#rule get_beagle:
#    output:
#        "beagle.21Apr21.304.jar"
#    log:
#        "logs/beagle/download.log"
#    resources: time_min=5, mem_mb=100, cpus=1
#    threads: 1
#    shell:
#        "wget http://faculty.washington.edu/browning/beagle/{output} > {log}"
#
#rule beagle_impute1:
#    input:
#        vcf="calls/all_merged.q20dp8.vcf.gz",
#        beagleexec="beagle.21Apr21.304.jar"
#    output:
#        "calls/all_merged.q20dp8.imputed.vcf.gz"
#    log:
#        "logs/beagle/beagle_imputefile2.log"
#    params:
#        prefix="calls/all_merged.q20dp8.imputed",
#        javamem="22"
#    resources: time_min=220, mem_mb=24000, cpus=16
#    threads: 16
#    conda:
#        "envs/bcftools.yaml"
#    shell:
#        "java -Xmx{params.javamem}g -jar {input.beagleexec} nthreads={threads} gt={input.vcf} out={params.prefix} 2> {log}"
#
#rule beagle_impute2:
#    input:
#        vcf="calls/all_merged.q20dp8saf0rpr1.vcf.gz",
#        beagleexec="beagle.21Apr21.304.jar"
#    output:
#        "calls/all_merged.q20dp8saf0rpr1.imputed.vcf.gz"
#    log:
#        "logs/beagle/beagle_imputefile2.log"
#    params:
#        prefix="calls/all_merged.q20dp8saf0rpr1.imputed",
#        javamem="22"
#    resources: time_min=820, mem_mb=24000, cpus=16
#    threads: 16
#    conda:
#        "envs/bcftools.yaml"
#    shell:
#        "java -Xmx{params.javamem}g -jar {input.beagleexec} nthreads={threads} gt={input.vcf} out={params.prefix} 2> {log}"
