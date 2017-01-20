<?php 

$variants_deva = $argv[1];
$serie = $argv[2];
$ngs_path = $argv[3];
$nm_file = $argv[4];

$patient_id = str_replace(".vcf", "", basename($variants_deva));


$variants_deva_snpeff = str_replace(".vcf", ".snpEff.vcf", $variants_deva);
$variants_deva_temp1 = str_replace(".vcf", ".temp1.vcf", $variants_deva);
$variants_deva_temp2 = str_replace(".vcf", ".temp2.vcf", $variants_deva);
$variants_deva_temp3 = str_replace(".vcf", ".temp3.vcf", $variants_deva);
$variants_deva_temp4 = str_replace(".vcf", ".temp4.vcf", $variants_deva);
$variants_deva_temp5 = str_replace(".vcf", ".temp5.vcf", $variants_deva);

system("java -Djava.io.tmpdir=$ngs_path/tmp -Xmx8g -jar $ngs_path/Programs/snpEff/SnpSift.jar annotate -v $ngs_path/Databases/snpEff/dbsnp_b147.vcf.gz $variants_deva > $variants_deva_temp1");

system("java -Djava.io.tmpdir=$ngs_path/tmp -Xmx8g -jar $ngs_path/Programs/snpEff/SnpSift.jar annotate -v $ngs_path/Databases/snpEff/refCosmic/Cosmic_v69_DeVa_Edition.vcf $variants_deva_temp1 > $variants_deva_temp2");

system("java -Djava.io.tmpdir=$ngs_path/tmp -Xmx8g -jar $ngs_path/Programs/snpEff/snpEff.jar eff -c $ngs_path/Programs/snpEff/snpEff.config -onlyTr $ngs_path/Genetics/$nm_file -v -hgvs -noLog -noStats -noMotif -noNextProt hg19 $variants_deva_temp2 > $variants_deva_temp3");

system("java -Djava.io.tmpdir=$ngs_path/tmp -Xmx8g -jar $ngs_path/Programs/snpEff/SnpSift.jar dbnsfp -f SIFT_pred,Polyphen2_HDIV_pred -v -db $ngs_path/Databases/snpEff/dbNSFP/dbNSFP/dbNSFP2.5.txt.gz $variants_deva_temp3 > $variants_deva_temp4");

system("java -Djava.io.tmpdir=/biopathdata/pipelineuser/NGS/tmp -Xmx8g -jar $ngs_path/Programs/snpEff/SnpSift.jar annotate -v $ngs_path/Databases/snpEff/ESP6500SI-V2-SSA137.GRCh38-liftover.vcf $variants_deva_temp4 > $variants_deva_temp5");
 
system("java -Djava.io.tmpdir=/biopathdata/pipelineuser/NGS/tmp -Xmx8g -jar $ngs_path/Programs/snpEff/SnpSift.jar annotate -v $ngs_path/Databases/snpEff/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz $variants_deva_temp5 > $variants_deva_snpeff");


unlink($variants_deva_temp1);
unlink($variants_deva_temp2);
unlink($variants_deva_temp3);
unlink($variants_deva_temp4);
unlink($variants_deva_temp5);

$dbh = new PDO("sqlite:$ngs_path/Genetics/RUNS/$serie/Analysis/DB/$serie.db");

$serie_rc = str_replace("-", "_", $serie);

$dbh->exec("CREATE TABLE IF NOT EXISTS variants_deva_$serie_rc (identifiant text, patient_id varchar(100), gene varchar(50), chr varchar(2), pos int, ref text, alt text, type varchar(10), refdepth int, altdepth int, allelicfreq double, exon varchar(100), codonchange text, aachange text, soft varchar(50), af_esp double, af_1000g double)");


if (file_exists($variants_deva_snpeff)) {
	$lines = file($variants_deva_snpeff);
	$count = count($lines);

	for($l = 0; $l < $count; $l++){
		$line = $lines[$l];
		if($line && ! preg_match("/^#/", $line)){
			$af_esp = $af_1000g = "";
			$line = trim($line);
			$elements = explode ("\t", $line);

			$chr = $elements[0];
			$pos = $elements[1];
			$ref = $elements[3];
			$alt = $elements[4];
			$type = (strlen($ref) == strlen($alt)) ? "SNV" : "INDEL";
			$info = $elements[9];
			$frags = explode (":", $info);
			$refalt = explode (",", $frags[1]);
			$refdepth = $refalt[0];
			$altdepth = $refalt[1];
			$allelicfreq = round($altdepth/($refdepth+$altdepth), 2);
			
			$more = $elements[7];

			$sgarf = explode (",", $more);

			if(count($sgarf) > 1){
				for($s = 0; $s < count($sgarf); $s++){
					if(preg_match("/^intron/", $sgarf[$s]) || preg_match("/UTR/", $sgarf[$s])){
						$more = $sgarf[$s];
					}
				}
			}

			$more = substr($more, strpos($more, ";EFF="));
			$stle = explode ("|", $more);
			$gene = $stle[5];
			$nm = $stle[8];
			$exon = $stle[9];
			$frags = explode("/", $stle[3]);
			$codonchange = (count($frags) == 2) ? $frags[1] : $frags[0];
			$aachange = (count($frags) == 2) ? $frags[0] : "";
		
			$soft = "deva";
			
			if(preg_match('/EA_AC=/', $more)){
				$af_esp = explode(",", explode("EA_AC=", $more)[1])[0]/(explode(",", explode("EA_AC=", $more)[1])[0]+explode(",", explode("EA_AC=", $more)[1])[1]);
			}
			
			if(preg_match('/EUR_AF=/', $more)){
				$af_1000g = substr(explode("EUR_AF=", $more)[1], 0, strpos(explode("EUR_AF=", $more)[1], ";"));
			}

			$info_deva = $patient_id."".$chr."".$pos."".$ref."".$alt;
			
			$dbh->exec("insert into variants_deva_$serie_rc values('$info_deva', '$patient_id', '$gene', '$chr', '$pos', '$ref', '$alt', '$type', '$refdepth', '$altdepth', '$allelicfreq', '$exon', '$codonchange', '$aachange', '$soft', '$af_esp', '$af_1000g')");
		}
	}
}

$dbh = null;
?>
