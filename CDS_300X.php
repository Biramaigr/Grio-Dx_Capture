<?php 

$bamfile = $argv[1];
$ngs_path = $argv[2];

$patient_id = str_replace(".bam", "", basename($bamfile));


echo "$patient_id;";

$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/BRCA1.bed");
$count = count($lines);
$bases_x300_brca1 = 0;
$bases_brca1 = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_brca1++;
		}
	}
	$bases_brca1 += abs($start-$stop)+1;
}
$bases_x300_percent_brca1 = 100*round($bases_x300_brca1/$bases_brca1, 2);
echo "$bases_x300_percent_brca1;";


$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/BRCA2.bed");
$count = count($lines);
$bases_x300_brca2 = 0;
$bases_brca2 = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_brca2++;
		}
	}
	$bases_brca2 += abs($start-$stop)+1;
}
$bases_x300_percent_brca2 = 100*round($bases_x300_brca2/$bases_brca2, 2);
echo "$bases_x300_percent_brca2;";


$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/ATR.bed");
$count = count($lines);
$bases_x300_atr = 0;
$bases_atr = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_atr++;
		}
	}
	$bases_atr += abs($start-$stop)+1;
}
$bases_x300_percent_atr = 100*round($bases_x300_atr/$bases_atr, 2);
echo "$bases_x300_percent_atr;";


$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/ATM.bed");
$count = count($lines);
$bases_x300_atm = 0;
$bases_atm = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_atm++;
		}
	}
	$bases_atm += abs($start-$stop)+1;
}
$bases_x300_percent_atm = 100*round($bases_x300_atm/$bases_atm, 2);
echo "$bases_x300_percent_atm;";


$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/NF1.bed");
$count = count($lines);
$bases_x300_nf1 = 0;
$bases_nf1 = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_nf1++;
		}
	}
	$bases_nf1+= abs($start-$stop)+1;
}
$bases_x300_percent_nf1 = 100*round($bases_x300_nf1/$bases_nf1, 2);
echo "$bases_x300_percent_nf1;";



$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/PTEN.bed");
$count = count($lines);
$bases_x300_pten = 0;
$bases_pten = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_pten++;
		}
	}
	$bases_pten+= abs($start-$stop)+1;
}
$bases_x300_percent_pten = 100*round($bases_x300_pten/$bases_pten, 2);
echo "$bases_x300_percent_pten;";


$lines = file("/biopathdata/pipelineuser/NGS/Genetics/CDS/TP53.bed");
$count = count($lines);
$bases_x300_tp53 = 0;
$bases_tp53 = 0;
for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	$line = trim($line);
	$elements = explode ("\t", $line);
	$chr = $elements[0];
	$start = ($elements[1]+10);
	$stop = ($elements[2]-10);
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	for($m = 0; $m < count($results); $m++){
		$dp = explode("\t", $results[$m])[2];
		if($dp >= 300){
			$bases_x300_tp53++;
		}
	}
	$bases_tp53+= abs($start-$stop)+1;
}
$bases_x300_percent_tp53 = 100*round($bases_x300_tp53/$bases_tp53, 2);
echo "$bases_x300_percent_tp53\n";
?>
