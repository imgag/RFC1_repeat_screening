# RFC1_repeat_screening

This is a tool for screening exome NGS data for the RFC1/CANVAS repeat.  
Because the repeat lies in an intron, it is typically not targeted by exome kits.  
However, we found that the locus is still covered by a few reads, at least for some enrichment kits e.g. SureSelect Human All Exon V6/7.

For WGS data or exome kits that target the repeat locus, more sophisticated methods like [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) should be used.

## Installation

The following dependencies need to be installed:

- php (for command line)
- tabix
- samtools

For Ubuntu 18.04 the installation of the dependencies can be done using this command:

	> sudo apt-get install php7.2-cli samtools tabix

After installing the dependencies, clone or download this repository and you are ready to go.


## Usage

After the installation you can call the tool like this:

	> php RFC1_repeat_screening.php [vcf_list] [gene_locus] [bam_list] [repeat_locus] > output.tsv

### Input

The tool needs the following input:

- *vcf_list:* A file containing the single-sample VCF.GZ files of the samples to screen (one per line).  
  Note: The VCF.GZ files have to contain the genotype *GT* as first entry of the *FORMAT* column!  
  Note: The VCF.GZ files have to be index!  
- *gene_locus:* The genomic coordinates of the RFC1 gene locus, i.e. `chr4:39287069-39370001` for GRCh37 and `chr4:39285449-39368381` for GRCH38.  
  Note: Remove the prefix *chr* if it is not contained in your reference genome. 
- *bam_list:* A file containing the BAM files of the samples to screen (one per line).  
  Note: The BAM files have to be index!  
  Note: The order of the VCF and BAM files has to be in sync, i.e. the first VCF file has to be from the same sample as the first BAM file, and so on.
- *repeat_locus:* The genomic coordinates of the RFC1 repeat locus, i.e. `chr4:39350034-39350137` for GRCh37 and `chr4:39348414-39348517` for GRCH38.  
  Note: Remove the prefix *chr* if it is not contained in your reference genome. 

### Output

The output is a tab-separated text file with one line per input samples.  
It has the following columns:

- *name:* Base file name of the BAM file.
- *variants:* Number of variants found in the VCF at the gene locus.
- *variants_hom_perc:* Percentage of above variants with homozygous genotype, i.e. homozygous alternate.
- *reads:* Number of reads found at the repeat locus.
- *soft_clipped_reads:* Number of soft-clipped reads found at the repeat locus.
- *AAGGG_repeats:* Number of AAGGG motif occurances in the soft-clipped reads parts.
- *other_repeats:* Other 5-mer motifs found in the soft-clipped reads parts.

## Filtering the output

A simple filtering strategy to select candidate samples for homozygous AAGGG repeat expansion would be: 

- variants >=1	
- variants_hom_perc >= 80	
- AAGGG_repeats >= 7	
- other_repeats empty


## How to cite

coming soon...

## License

This tool is provided under [MIT license](LICENSE).
