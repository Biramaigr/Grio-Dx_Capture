<?php

$serie = $argv[1];
$ngs_path = $argv[2];
$mode = $argv[3];
$bedfile = $argv[4];

echo "patient_id;nb_reads_fw;nb_reads_q30_fw;%read_q30_fw;nb_reads_rv;nb_reads_q30_rv;%read_q30_rv;read_length_mean_fw;read_quality_mean_fw;read_length_mean_rv;read_quality_mean_rv;nb_mapped_reads_fw;%mapped_reads_fw;nb_mapped_reads_remove_dup_fw;%mapped_reads_remove_dup_fw;nb_mapped_reads_rv;%mapped_reads_rv;nb_mapped_reads_remove_dup_rv;%mapped_reads_remove_dup_rv;cds_percent_100;cds_percent_300;cds_percent_ard_100;cds_percent_ard_300\n";

if ($handle = opendir("$ngs_path/Genetics/RUNS/$serie/Analysis/QCAs_tech")) {
    while (false !== ($entry = readdir($handle))) {
	if (preg_match("/.QCA_tech.csv$/", $entry)) {
      		$patient_id = str_replace(".QCA_tech.csv", "", $entry);
		$patient_id1 = substr($patient_id, 0, strlen($patient_id)-4);
		
		system("php $ngs_path/Scripts/Genetics_Capture/$mode/QC_Library.php $ngs_path/Genetics/RUNS/$serie/fastQ/".$patient_id1."_R1_001.fastq $ngs_path/Genetics/RUNS/$serie/fastQ/".$patient_id1."_R2_001.fastq $ngs_path/Genetics/RUNS/$serie/Analysis/Bams_RAW/".$patient_id.".rdup.bam $ngs_path/Genetics/RUNS/$serie/Analysis/Bams_RAW/".$patient_id.".nrr.bam $bedfile $ngs_path");
	} 
    }
    closedir($handle);
}

?>
