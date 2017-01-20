<?php 

$bamfile = $argv[1];
$bedfile = $argv[2];
$serie = $argv[3];
$ngs_path = $argv[4];
$dp_reseq = $argv[5];

include("$ngs_path/Programs/pChart/class/pData.class.php");
include("$ngs_path/Programs/pChart/class/pDraw.class.php");
include("$ngs_path/Programs/pChart/class/pImage.class.php");
include("$ngs_path/Programs/pChart/class/pPie.class.php");
include("$ngs_path/Programs/pChart/class/pIndicator.class.php");

$patient_id = str_replace(".bam", "", basename($bamfile));

$lines = file($bedfile);
$count = count($lines);

$fpn = fopen("$ngs_path/Genetics/RUNS/$serie/Analysis/QCAs_DP_Normalize/".$patient_id.".QCA_DP_Normalize.csv", 'w');

fwrite($fpn, "patient_id;amplicon;chr;start;stop;dp_mean;status;exon;dm_min;dp_max\n");

$array_x = array();
$array_y = array();

for($l = 0; $l < $count; $l++){
	$results = null;
	$line = $lines[$l];
	#echo "$line";
	$line = trim($line);
	$elements = explode ("\t", $line);

	$amplicon = $elements[0];
	$chr = $elements[1];
	$start = $elements[2];
	$stop = $elements[3];
	$exon = $elements[4];
	
	exec("$ngs_path/Programs/samtools-1.2/samtools depth -r $chr:$start-$stop $bamfile", $results);
	
	$array_dp = array();
	$summ = 0;
	for($m = 0; $m < count($results); $m++){
		$summ += explode("\t", $results[$m])[2];
		array_push($array_dp, explode("\t", $results[$m])[2]);
	}
	
	$results_size = (count($results) == 0) ? 1 : count($results);

	$dp_mean = round($summ/$results_size);
	
	if($dp_mean >= $dp_reseq){
		$status = "PASS";
	}
	else{
		$status = "FAIL";
	}
	
	fwrite($fpn, "$patient_id;$amplicon;$chr;$start;$stop;$dp_mean;$status;$exon;".min($array_dp).";".max($array_dp)."\n");
}

fclose($fpn);


?>
