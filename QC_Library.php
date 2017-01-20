<?php

$fastq_r1 = $argv[1];
$fastq_r2 = $argv[2];
$bamfile = $argv[3];
$bamfile_raw = $argv[4];
$bedfile = $argv[5];
$ngs_path = $argv[6];

$patient_id = str_replace(".rdup.bam", "", basename($bamfile));

$lines = file($fastq_r1);
$count = count($lines);

$nb_reads_fw = $count/4;
$nb_reads_q30_fw = 0;
$read_quality_global_tot_fw = 0;

$read_length_tot_fw = 0;

$array1 = array();

for($l = 0; $l < $count; $l++){
	$read_id = trim($lines[$l]);
	$read_seq = trim($lines[$l+1]);
	$read_qual = trim($lines[$l+3]);

	$read_length = (strlen($read_seq) > 0) ? strlen($read_seq) : 1;

	array_push($array1, $read_length);

	$read_length_tot_fw += $read_length;

	$read_quality_tot = 0;
	for($q = 0; $q < $read_length; $q++){
		$char_current = substr($read_qual, $q, 1);
		$qual_current = ord($char_current)-33;
		$read_quality_tot += $qual_current;
	}

	$read_quality_mean = round($read_quality_tot/$read_length);

	$read_quality_global_tot_fw += $read_quality_mean;

	if($read_quality_mean >= 30){
		$nb_reads_q30_fw++;
	}

	$l = $l + 3;
}

$read_length_mean_fw = round(4*$read_length_tot_fw/$count);

$min_fw = min($array1);
$max_fw = max($array1);
$read_length_distrib_fw = "$read_length_mean_fw [$min_fw<->$max_fw]";
$read_quality_global_mean_fw = round(4*$read_quality_global_tot_fw/$count);


$lines2 = file($fastq_r2);
$count2 = count($lines2);

$nb_reads_rv = $count2/4;
$nb_reads_q30_rv = 0;
$read_quality_global_tot_rv = 0;

$read_length_tot_rv = 0;

$array2 = array();

for($l = 0; $l < $count2; $l++){
	$read_id2 = trim($lines2[$l]);
	$read_seq2 = trim($lines2[$l+1]);
	$read_qual2 = trim($lines2[$l+3]);

	$read_length2 = (strlen($read_seq2) > 0) ? strlen($read_seq2) : 1;

	array_push($array2, $read_length2);

	$read_length_tot_rv += $read_length2;

	$read_quality_tot2 = 0;
	for($q = 0; $q < $read_length2; $q++){
		$char_current2 = substr($read_qual2, $q, 1);
		$qual_current2 = ord($char_current2)-33;
		$read_quality_tot2 += $qual_current2;
	}

	$read_quality_mean2 = round($read_quality_tot2/$read_length2);

	$read_quality_global_tot_rv += $read_quality_mean2;

	if($read_quality_mean2 >= 30){
		$nb_reads_q30_rv++;
	}

	$l = $l + 3;
}

$read_length_mean_rv = round(4*$read_length_tot_rv/$count2);

$min_rv = min($array2);
$max_rv = max($array2);
$read_length_distrib_rv = "$read_length_mean_rv [$min_rv<->$max_rv]";
$read_quality_global_mean_rv = round(4*$read_quality_global_tot_rv/$count2);


exec("$ngs_path/Programs/samtools-1.2/samtools view -F 16 $bamfile_raw | wc -l", $results_raw1);
exec("$ngs_path/Programs/samtools-1.2/samtools view -f 16 $bamfile_raw | wc -l", $results_raw2);

$nb_mapped_reads_raw_fw = $results_raw1[0];
$nb_mapped_reads_raw_rv = $results_raw2[0];

$mapped_reads_percent_raw_fw = round(100*$nb_mapped_reads_raw_fw/$nb_reads_fw);
$mapped_reads_percent_raw_rv = round(100*$nb_mapped_reads_raw_rv/$nb_reads_rv);


exec("$ngs_path/Programs/samtools-1.2/samtools view -F 16 $bamfile | wc -l", $results1);
exec("$ngs_path/Programs/samtools-1.2/samtools view -f 16 $bamfile | wc -l", $results2);

$nb_mapped_reads_fw = $results1[0];
$nb_mapped_reads_rv = $results2[0];

$mapped_reads_percent_fw = round(100*$nb_mapped_reads_fw/$nb_reads_fw);
$mapped_reads_percent_rv = round(100*$nb_mapped_reads_rv/$nb_reads_rv);

$cds_percent_100 = getCDSPercent($bamfile_raw, $bedfile, 100, $ngs_path);
$cds_percent_300 = getCDSPercent($bamfile_raw, $bedfile, 300, $ngs_path);
$cds_percent_ard_100 = getCDSPercent($bamfile, $bedfile, 100, $ngs_path);
$cds_percent_ard_300 = getCDSPercent($bamfile, $bedfile, 300, $ngs_path);

echo "$patient_id;$nb_reads_fw;$nb_reads_q30_fw;".round(100*$nb_reads_q30_fw/$nb_reads_fw).";$nb_reads_rv;$nb_reads_q30_rv;".round(100*$nb_reads_q30_rv/$nb_reads_rv).";$read_length_distrib_fw;$read_quality_global_mean_fw;$read_length_distrib_rv;$read_quality_global_mean_rv;$nb_mapped_reads_raw_fw;$mapped_reads_percent_raw_fw;$nb_mapped_reads_fw;$mapped_reads_percent_fw;$nb_mapped_reads_raw_rv;$mapped_reads_percent_raw_rv;$nb_mapped_reads_rv;$mapped_reads_percent_rv;$cds_percent_100;$cds_percent_300;$cds_percent_ard_100;$cds_percent_ard_300\n";


function getCDSPercent($bamfile, $bedfile, $dp_min, $ngs_path){

	$lines = file("$bedfile");
	$count = count($lines);
	$bases_xDP = 0;
	$bases = 0;
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
			if($dp >= $dp_min){
				$bases_xDP++;
			}
		}
		$bases += abs($start-$stop)+1;
	}
	$bases_xDP_percent = 100*round($bases_xDP/$bases, 2);

	return $bases_xDP_percent;

}

?>
