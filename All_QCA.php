<?php

$bamdir = $argv[1];
$serie = $argv[2];
$ngs_path = $argv[3];
$mode = $argv[4];
$bed_name = $argv[5];
$dp_reseq = $argv[6];
$dp_conta = $argv[7];

if ($handle = opendir($bamdir)) {
    while (false !== ($entry = readdir($handle))) {
        if (preg_match("/.bam$/", $entry)) {
			system("php $ngs_path/Scripts/Genetics_Capture/$mode/QCA.php $bamdir/$entry $ngs_path/Genetics/$bed_name $serie $ngs_path $dp_reseq");

			$bedfile_dp_normalize = str_replace("_name.bed", "_DP_Normalize.bed", $bed_name);
			
			system("php $ngs_path/Scripts/Genetics_Capture/$mode/QCA_DP_Normalize.php $bamdir/$entry $ngs_path/Genetics/$bedfile_dp_normalize $serie $ngs_path $dp_reseq");

			system("php $ngs_path/Scripts/Genetics_Capture/$mode/QCA_tech.php $bamdir/$entry $ngs_path/Genetics/$bed_name $serie $ngs_path $dp_reseq $dp_conta");
        }
    }
    closedir($handle);
}

?>
