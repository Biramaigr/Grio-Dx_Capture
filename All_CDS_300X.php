<?php

$serie = $argv[1];
$ngs_path = $argv[2];
$mode = $argv[3];

echo "patient_id;brca1;brca2;atr;atm;nf1;pten;tp53\n";

if ($handle = opendir("$ngs_path/Genetics/RUNS/$serie/Analysis/QCAs_tech")) {
    while (false !== ($entry = readdir($handle))) {
	if (preg_match("/.QCA_tech.csv$/", $entry)) {
      		$patient_id = str_replace(".QCA_tech.csv", "", $entry);
		$patient_id1 = substr($patient_id, 0, strlen($patient_id)-4);
		
		system(" php $ngs_path/Scripts/Genetics_Capture/$mode/CDS_300X.php $ngs_path/Genetics/RUNS/$serie/Analysis/Bams/".$patient_id.".bam $ngs_path");
	} 
    }
    closedir($handle);
}

?>
