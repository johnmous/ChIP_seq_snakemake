## Author: I. Moustakas, i.m*@lumc.nl
## Title: Analysis to call differential peaks for 102 ratChIP-seq
## IMPORTANT NOTICE: We use the old gff file as input to HTseq. All steps up to there are thus irrelevant for the output

import yaml

with open("/path/to/sample_sheet_controls.yml", "r") as config: 
    try:
        configDict = yaml.load(config)
    except yaml.YAMLError as exc:
        print(exc)

narrowPeakList = []
for key, value in configDict['samples'].items():
    if 'control' in value.keys():
        control = value["control"]
        sample = key
        narrowPeakList.append("/{0}/macs2/{0}_VS_{1}/{0}_VS_{1}_peaks.narrowPeak".format(sample, control))



SAMPLES = narrowPeakList

INPUT_PATH="/path/to/project-102-ratchipseq/analysis/run4_carp_Rn_6/samples"
OUTPUT_PATH="/path/to/project-102-ratchipseq/analysis/run4_carp_Rn_6/diff_analysis_dedup"

## A function to grep for a substring in a list of strings. Returns the list of strings containing substring
def grepList(listFiles, substr):
    listOfFiles = list(filter(lambda x: substr in x, listFiles))
    return(listOfFiles)

## A function to get the sample IDs from the narrowPeak list of files after filtering for Transcription Factor
def getSampleID(listNarrowPeak, TF):
    sampleIDs=[]
    for file in listNarrowPeak:
        sampleIDs.append(file.split("/")[1])
    sampleIDsPerTF = grepList(sampleIDs, TF)
    return sampleIDsPerTF

rule all: 
    input: 
        expand(OUTPUT_PATH + "/{TF}_tablesMerged.count", TF=["GR", "pCREB"])
#         expand(OUTPUT_PATH + "/{TF}_counts/{sampleID}.count", TF=["GR", "pCREB"])

narrowPeakFiles = [INPUT_PATH + s for s in SAMPLES]
## Concatente the narrow peak files 
## Lambda function as input to get all samples for one TF
## mergeBed: collapse the peak names (column 4), sum the peak scores (column 5)
rule concatenate:
    input: 
        lambda wildcards: grepList(narrowPeakFiles, wildcards.TF)
    output: 
        OUTPUT_PATH + "/{TF}_merged.bed"
    params:
        sthreads=1,
        mem="2G",
        time="1:0:0"
    conda:
       "/path/to/chipseq.yml"    
    shell:
        """
        cat {input} | sort -k1,1 -k2,2n | mergeBed -i stdin -c 4,5 -o collapse,sum  > {output}
        """

## convert bed to gff using an R script
rule convertBed2Gff:
    input:
        OUTPUT_PATH + "/{TF}_merged.bed"
    output:
         OUTPUT_PATH + "/{TF}_merged.gff"
    params:
        sthreads=1,
        mem="4G",
        time="1:0:0"
    conda:
       "/path/to/chipseq.yml"    
    shell:
        """
        Rscript /path/to/get_gff.R {input} {output}
        """

## Replace asterisc with dot for column 7 of the GFF file
rule fixGFF:
    input:
        OUTPUT_PATH + "/{TF}_merged.gff"
    output:
        OUTPUT_PATH + "/{TF}_mergedFixAsterisc.gff"
    params:
        sthreads=1,
        mem="4G",
        time="1:0:0"
    shell:
        """
        cat {input} | awk  '{{if ($7 == "*") {{sub($7, "."); print }} else print }}' > {output}
        """


## Shorten peak IDs
rule shortenPeakIDs:
    input:
        OUTPUT_PATH + "/{TF}_mergedFixAsterisc.gff"
    output:
        OUTPUT_PATH + "/{TF}_merged_IDsShort.gff"
    params:
        sthreads=1,
        mem="1G",
        time="1:0:0"
    shell: 
        """
        cat {input} | sed 's/_VS_S_[0-9]\+_Input_[0-9]\+_peak//g' > {output}
        """

## picard tools remove duplicates
rule picard:
    input: 
        INPUT_PATH + "/{sampleID}/{sampleID}.filter.bam"
    output: 
        INPUT_PATH + "/{sampleID}/{sampleID}.filter.dedupl.bam"
    params: 
        metrics=INPUT_PATH + "/{sampleID}/{sampleID}.marked_dup_metrics.txt",
        sthreads=2,
        mem="8G",
        time="6:0:0"
    conda:
        "/path/to/chipseq.yml"
    shell:
        """
        //path/to/miniconda3/envs/chipseq/share/picard-2.23.9-0/picard -Xmx6656M MarkDuplicates \
        REMOVE_DUPLICATES=true \
        I={input} \
        O={output} \
        M={params.metrics}
        """


## Count the number of reads mapping on each peak. 
## htseq-counts settings: 
## -s no, non stranded RNA experiment
## -f bam, input file is bam format
## -t sequence_feature, feature type to count on (in our case all peaks are named sequence_feature)
## -m intersection-strict, mode to handle reads overlapping more than one feature is strict
rule countPeaks:
    input: 
        bam=INPUT_PATH + "/{sampleID}/{sampleID}.filter.dedupl.bam",
        gff="/path/to/project-102-ratchipseq/analysis/run4_carp_Rn_6/diffAnalysis/{TF}_merged_IDsShort.gff"
        #OUTPUT_PATH + "/{TF}_merged_IDsShort.gff"
    output: 
        OUTPUT_PATH + "/{TF}_counts/{sampleID}.count"
    params: 
        sthreads=1,
        mem="12G",
        time="3:0:0"
    conda:
        "/path/to/chipseq.yml"
    shell:
        """
        htseq-count -s no -m intersection-strict -t sequence_feature -f bam {input.bam} {input.gff} > {output}
        """ 


## Merge the count tables
## Use a lambda function to loop over the TF and get the appropriate sampleIDs per TF
rule mergeTables:
    input: 
        lambda wildcards: expand("{outputPath}/{TF}_counts/{sampleID}.count" , sampleID = getSampleID(SAMPLES, wildcards.TF), outputPath=OUTPUT_PATH, TF = wildcards.TF)
    output:
        OUTPUT_PATH + "/{TF}_tablesMerged.count"
    params: 
        sthreads=1,
        mem="10G",
        time="3:0:0"
    shell:
        """
        collect-columns {output} {input}
        """


