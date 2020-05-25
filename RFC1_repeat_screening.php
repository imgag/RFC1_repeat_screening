<?php

################################## check input ##################################
if ($argc!=5)
{
	fwrite(STDERR, "Usage: RFC1_repeat_screening.php [vcf_list] [gene_locus] [bam_list] [repeat_locus]\n");
	exit(1);
}

$vcfs = file($argv[1]);
if ($vcfs===FALSE)
{
	fwrite(STDERR, "Error: VCF list file '".$argv[1]."' is not readable!\n");
	exit(1);
}
$vcfs = array_map("trim", $vcfs);
$vcfs = array_filter($vcfs);

$gene_locus = $argv[2];

$bams = file($argv[3]);
if ($bams===FALSE)
{
	fwrite(STDERR, "Error: BAM list file '".$argv[3]."' is not readable!\n");
	exit(1);
}
$bams = array_map("trim", $bams);
$bams = array_filter($bams);

if (count($vcfs)!=count($bams))
{
	fwrite(STDERR, "Error: VCF and BAM count in input files differs!\n");
	exit(1);
}

$repeat_locus = $argv[4];


################################## auxilary functions ##################################
function correct($mer)
{
	$shifts = array();
	for ($i=0; $i<strlen($mer); ++$i)
	{
		$shifts[] = substr($mer.$mer, $i, 5);
	}
	sort($shifts);	
	return $shifts[0];
}

function exec2($command)
{
    //start processing
	$proc = proc_open($command, array(1 => array('pipe','w'), 2 => array('pipe','w')), $pipes);
	
	//get stdout, stderr and exit code
    $stdout = stream_get_contents($pipes[1]);
    fclose($pipes[1]);
    $stderr = stream_get_contents($pipes[2]);
    fclose($pipes[2]);
    $exit = proc_close($proc);
	
	//abort if requested
	if ($exit!=0)
	{
		trigger_error("Error while executing command: '$command'\nCODE: $exit\nSTDOUT: ".$stdout."\nSTDERR: ".$stderr."\n", E_USER_ERROR);
	}
	
	//return output
	return array(explode("\n", rtrim($stdout)), explode("\n", rtrim($stderr)), $exit);
}

################################## screening ##################################
print "#name\tvariants\tvariants_hom_perc\treads\tsoft_clipped_reads\tAAGGG_repeats\tother_repeats\n";
for ($i=0; $i<count($vcfs); ++$i)
{
	$vcf = $vcfs[$i];
	if (!file_exists($vcf))
	{
		fwrite(STDERR, "Error: VCF file '$vcf' is not readable!\n");
		exit(1);
	}
	
	$bam = $bams[$i];
	if (!file_exists($vcf))
	{
		fwrite(STDERR, "Error: VCF file '$vcf' is not readable!\n");
		exit(1);
	}

	//variant genotypes
	list($vars) = exec2("tabix {$vcf} {$gene_locus}");
	$var_count = 0;
	$hom_count = 0;
	foreach($vars as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		list($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample) = explode("\t", $line);
		++$var_count;
		
		list($gt) = explode(":", $sample);
		if ($gt=="1/1" || $gt=="1|1") ++$hom_count;
	}
		
	//determine reads
	list($reads) = exec2("samtools view {$bam} {$repeat_locus}");
	
	//determine soft-clipped repeat
	$dp = 0;
	$soft_clipped = 0;
	$rep_AAGGG = 0;
	$five_mers = array();
	foreach($reads as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		
		++$dp;
		
		list($id, $flags, $chr, $pos, $mq, $cigar, , , , $bases) = explode("\t", $line);
		if (strpos($cigar, "S")!==FALSE)
		{
			$matches = array();
			preg_match_all("/[0-9]+[A-Z]/", $cigar, $matches);
			foreach($matches[0] as $match)
			{
				$len = substr($match, 0, strlen($match)-1); 
				$op = substr($match, strlen($match)-1 , 1);
				$seq = substr($bases, 0, $len);
				$bases = substr($bases, $len);
				if ($op=="S")
				{
					$rep_AAGGG += substr_count($seq, "AAGGG");
					for($j=0; $j<strlen($seq)-4; ++$j)
					{
						@$five_mers[correct(substr($seq, $j, 5))] += 1;
					}
				}
			
			}
			 ++$soft_clipped;
		}
	}
	$hom_frac = $var_count==0 ? "n/a" : number_format(100.0*$hom_count/$var_count, 2);
	$rep_other = array();
	arsort($five_mers);
	foreach($five_mers as $mer => $count)
	{
		if ($mer=="AAAAG") continue; //WT
		if ($mer=="AAGGG") continue; //pahogenic repeat
		$count = (int)($count/5);
		if ($count > 20) $rep_other[] = "$count*$mer";	
	}
	
	$name = basename($bam, ".bam");
	print "{$name}\t{$var_count}\t{$hom_frac}\t{$dp}\t{$soft_clipped}\t{$rep_AAGGG}\t".implode(", ", $rep_other)."\n";
}

?>