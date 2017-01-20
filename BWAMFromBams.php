<?php

$serie = $argv[1];
$ngs_path = $argv[2];
$mode = $argv[3];
$bams_group = $argv[4];
$bed = $argv[5];
$depth_variation_threshold = $argv[6];
$af_threshold = $argv[7];
$dp_threshold = $argv[8];

$files = getFiles("$bams_group");

sort($files);

for($f = 0; $f < count($files); $f++){

	$bamfile = "$bams_group/".$files[$f];
	$out = preg_replace("/.bam/", "",$files[$f]);

	system("samtools index $bamfile");

	//system("java -Djava.io.tmpdir=$ngs_path/tmp -Xmx8g -jar $ngs_path/Programs/GATK/v-3.4-46/GenomeAnalysisTK.jar -T HaplotypeCaller -R $ngs_path/Databases/hg19_chr/hg19.fa -I $bamfile --genotyping_mode DISCOVERY -stand_emit_conf 0 -stand_call_conf 0 -L $ngs_path/Genetics/$bed -o $ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_deva/$out.vcf -rf BadCigar");
	//unlink("$ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_deva/$out.vcf.idx");

	system("java -Djava.io.tmpdir=$ngs_path/tmp -Xmx8g -jar $ngs_path/Programs/MutaCaller-1.6/MutaCaller-1.6-rt.jar -LIBS_PATH=$ngs_path/Programs/MutaCaller-1.6 -MATE_TYPE=PAIRED -DPMIN=$dp_threshold -AFMIN=$af_threshold -PQUALMIN=10 -CLUSTER_WINDOW=20 -CLUSTER_NUMBER=4 -REFERENCE=$ngs_path/Databases/hg19_chr/hg19.fa -BEDFILE=$ngs_path/Genetics/$bed -BIG_INDELS=NO -BIG_INDELS_OPTIONS=$depth_variation_threshold,10,10,0.15,YES -INPUT=$bamfile -OUTPUT=$ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_mutacaller/$out.ns.vcf 2>>$ngs_path/Genetics/RUNS/$serie/Analysis/Logs/$out.log");
	system("cat $ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_mutacaller/$out.ns.vcf | $ngs_path/Programs/vcftools_0.1.12b/bin/vcf-sort > $ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_mutacaller/$out.vcf");
	unlink("$ngs_path/Genetics/RUNS/$serie/Analysis/Resultats_mutacaller/$out.ns.vcf");
	
	mkdir("$ngs_path/Genetics/RUNS/$serie/Analysis/Synchro/$out");
}

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

?>
