<?php

$serie = $argv[1];
$bwa_algo = $argv[2];
$is_pp = $argv[3];
$ngs_path = $argv[4];
$mode = $argv[5];
$length_min = $argv[6];
$is_indel_realign = $argv[7];
$db_extension = $argv[8];
$bed = $argv[9];
$bed_name = $argv[10];
$nm_file = $argv[11];
$dp_reseq = $argv[12];
$dp_conta = $argv[13];
$depth_variation_threshold = $argv[14];
$af_threshold = $argv[15];
$dp_threshold = $argv[16];

if(strtoupper($mode) == "PRODUCTION"){
	$mode = "Production";
}
else{
	$mode = "Development";
}

$run = $serie;
$serie = str_replace("-", "_", $serie);

mkdir("$ngs_path/Genetics/RUNS/$serie");
mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis");
mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Bams");
//mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group1");
//mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group2");
//mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group3");
//mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group4");
mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_deva");
mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_mutacaller");
mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Logs");
mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Synchro");

system("cp /biopathnas/NGS_Génétique/RUNS/$run/*.bam $ngs_path/Genetics/RUNS/$serie/Analysis/Bams/");

//$size = GroupPartition("$ngs_path/Genetics/RUNS/$serie/Analysis/Bams", "life");

system("php $ngs_path/Scripts/Genetics_Capture/$mode/BWAMFromBams.php $serie $ngs_path $mode $ngs_path/Genetics/RUNS/$serie/Analysis/Bams $bed $depth_variation_threshold $af_threshold $dp_threshold 1> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/output.txt 2> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/error.txt");

/*system("php $ngs_path/Scripts/Genetics_Capture/$mode/BWAMFromBams.php $serie $ngs_path $mode $ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group2 $bed $depth_variation_threshold $af_threshold $dp_threshold 1> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/output.txt 2> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/error.txt&");
system("php $ngs_path/Scripts/Genetics_Capture/$mode/BWAMFromBams.php $serie $ngs_path $mode $ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group3 $bed $depth_variation_threshold $af_threshold $dp_threshold 1> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/output.txt 2> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/error.txt&");
system("php $ngs_path/Scripts/Genetics_Capture/$mode/BWAMFromBams.php $serie $ngs_path $mode $ngs_path/Genetics/RUNS/$serie/Analysis/Bams_group3 $bed $depth_variation_threshold $af_threshold $dp_threshold 1> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/output.txt 2> $ngs_path/Genetics/RUNS/$serie/Analysis/Logs/error.txt&");

while(1 > 0){
	if(alreadydone("$ngs_path/Genetics/RUNS/$serie/Analysis/Synchro") == $size){
		break;
	}
}
*/
#IGR DATABASE COMMIT
system("cp $ngs_path/Databases/igr_apache.db $ngs_path/Databases/igr.db");
system("php $ngs_path/Scripts/Genetics_Capture/$mode/FinalStepFromBams.php $serie $ngs_path $mode $bed_name $dp_reseq $dp_conta $nm_file $bed");

system("php $ngs_path/Scripts/Genetics_Capture/$mode/DB_Storage.php $ngs_path/Genetics/RUNS/$serie/Analysis/Synthesis $ngs_path/Genetics/RUNS/$serie/Analysis/QCAs $ngs_path/Genetics/RUNS/$serie/Analysis/QC_Libraries.csv $ngs_path/Genetics/RUNS/$serie/Analysis/CDS_300X.csv $serie $run $ngs_path $mode $db_extension");

system("cat $ngs_path/Genetics/RUNS/$serie/Analysis/Synthesis/* > $ngs_path/Genetics/RUNS/$serie/Analysis/Synthesis/AllPatients.csv");


system("rm -rf $ngs_path/Genetics/RUNS/$serie/Analysis/Bams_RAW");
system("rm -rf $ngs_path/Genetics/RUNS/$serie/Analysis/Synchro");
system("rm -rf $ngs_path/Genetics/RUNS/$serie/fastQ_group1");
system("rm -rf $ngs_path/Genetics/RUNS/$serie/fastQ_group2");
system("rm -rf $ngs_path/Genetics/RUNS/$serie/fastQ_group3");
system("rm -rf $ngs_path/Genetics/RUNS/$serie/fastQ_group4");
system("rm -rf $ngs_path/Genetics/RUNS/$serie/fastQ");

mkdir("/biopathnas/NGS_Génétique/Analyses/$serie");
system("cp -r $ngs_path/Genetics/RUNS/$serie/Analysis/* /biopathnas/NGS_Génétique/Analyses/$serie/");


function getFiles($dir){
	$listfiles = scandir($dir);
	$tab = array();
	for($f = 0; $f < count($listfiles); $f++){
		$entry = $listfiles[$f];
		
		if (!preg_match("/^\./", $entry)) {
			array_push($tab, $entry);	
		}
	}
	
	sort($tab);
	return $tab;
}


function alreadydone($dir1) {
	$elts = array();
	if ($handle = opendir($dir1)) {
    	while (false !== ($entry = readdir($handle))) {
    		if ($entry != "." && $entry != ".."){
        		array_push($elts, $entry);
        	}
    	}
    	closedir($handle);
	}
	
	$size = count($elts);
	
	return $size;
}

function GroupPartition($folder, $sequencer){
	
	$elts = getFiles($folder);
	$increment = 0;
	$nbtot = count($elts);

	for($f = 0; $f < count($elts); $f++){
		$entry = $elts[$f];
	
		$increment++;
		
		if($nbtot < 8){
			if($sequencer == "illumina"){
				if($nbtot == 2){
					system("cp $folder/$entry ".$folder."_group1/");
				}
				else if($nbtot == 4){
					if($increment <= 2){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group2/");
					}
				}
				else if($nbtot == 6){
					if($increment <= 2){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment <= 4){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group3/");
					}
				}
			}
			else{
				if($nbtot == 1){
					system("cp $folder/$entry ".$folder."_group1/");
				}
				else if($nbtot == 2){
					if($increment == 1){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group2/");
					}
				}
				else if($nbtot == 3){
					if($increment == 1){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment == 2){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group3/");
					}
				}
				else if($nbtot == 4){
					if($increment == 1){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment == 2){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else if($increment == 3){
						system("cp $folder/$entry ".$folder."_group3/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group4/");
					}
				}
				else if($nbtot == 5){
					if($increment == 1){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment == 2){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else if($increment == 3){
						system("cp $folder/$entry ".$folder."_group3/");
					}
					else if($increment == 4){
						system("cp $folder/$entry ".$folder."_group4/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group1/");
					}
				}
				else if($nbtot == 6){
					if($increment == 1){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment == 2){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else if($increment == 3){
						system("cp $folder/$entry ".$folder."_group3/");
					}
					else if($increment == 4){
						system("cp $folder/$entry ".$folder."_group4/");
					}
					else if($increment == 5){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group2/");
					}
				}
				else if($nbtot == 7){
					if($increment == 1){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment == 2){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else if($increment == 3){
						system("cp $folder/$entry ".$folder."_group3/");
					}
					else if($increment == 4){
						system("cp $folder/$entry ".$folder."_group4/");
					}
					else if($increment == 5){
						system("cp $folder/$entry ".$folder."_group1/");
					}
					else if($increment == 6){
						system("cp $folder/$entry ".$folder."_group2/");
					}
					else{
						system("cp $folder/$entry ".$folder."_group3/");
					}
				}
			}
		}
		else{
			
			if($increment <= (intval($nbtot/4) + (intval($nbtot/4))%2)){
				system("cp $folder/$entry ".$folder."_group1/");
			}
			else if($increment <= intval($nbtot/4) + (intval($nbtot/4))%2 + intval($nbtot/4) + (intval($nbtot/4))%2){
				system("cp $folder/$entry ".$folder."_group2/");
			}
			else if($increment <= intval($nbtot/4) + (intval($nbtot/4))%2 + intval($nbtot/4) + (intval($nbtot/4))%2 + intval($nbtot/4) + (intval($nbtot/4))%2){
				system("cp $folder/$entry ".$folder."_group3/");
			}
			else{
				system("cp $folder/$entry ".$folder."_group4/");
			}
		}
	}
	
	return $nbtot;

}

?>
