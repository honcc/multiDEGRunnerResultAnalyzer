#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Math::Random;
use Math::Round;
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use Data::Dumper::Names;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to run GOTermFinder, GSEA, HyperFUNC, GOSlimmer and RamiGO, on the differential expression results specified in geneDEResultPathFile;
#
#	Input
#		--geneDEResultPathFile=	file path; path contains the Differential Expression results generated from multiCuffDiffDeSeqRunner;
#		--geneInfoPath=			file path; path of the whole genome gene ID and description; one line one gene format;
#		--outDir=				directory for output;
#
#		v0.1
#			[20/11/2013 14:45] renamed from multiCuffDiffDeSeqRunnerResultAnanlyzer_v0.3
#			[20/11/2013 14:53] removed runGOTermFinderInAllComparison
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-11-20 14:46]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GeneOntology/multiDEGRunnerResultAnalyzer/v0.1/multiDEGRunnerResultAnalyzer_v0.1.pl --geneInfoPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/GFFSeqExtractor/log/gene.info.tsv --geneDEResultPathFile=/Volumes/A_MPro2TB/NGS/analyses/EHI_heatShockDiffExp/All2Rep/HM1HeatShockAll2Rep.allDEGList.txt
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GeneOntology/multiDEGRunnerResultAnalyzer/v0.1/multiDEGRunnerResultAnalyzer_v0.1.pl
#	--geneInfoPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/GFFSeqExtractor/log/gene.info.tsv
#	--geneDEResultPathFile=/Volumes/A_MPro2TB/NGS/analyses/EHI_heatShockDiffExp/All2Rep/HM1HeatShockAll2Rep.allDEGList.txt
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $ARGVStr = join "\n", (&currentTime(), abs_path($0), @ARGV);#->280
my $globalScriptDirPath = dirname(rel2abs($0));
open DEBUGLOG, ">", "$globalScriptDirPath/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|818, readParameters|1140
#	secondaryDependOnSub: currentTime|280
#
#<section ID="startingTasks" num="0">
#----------Read parameters ----------#
&printCMDLogOrFinishMessage("CMDLog");#->818
my ($geneInfoPath, $geneDEResultPathFile, $outDir) = &readParameters();#->1140
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="1">
my @mkDirAry;

my $GSEADir = "$outDir/GSEA/"; push @mkDirAry, $GSEADir;
my $FUNCHyperDir = "$outDir/FUNCHyper/"; push @mkDirAry, $FUNCHyperDir;
my $dirForGSEAMasterHtml = $GSEADir;
my $dirForFUNCMasterHtml = $FUNCHyperDir;
my $geneBasedSummaryOutdir = "$outDir/geneBasedSummary/"; push @mkDirAry, $geneBasedSummaryOutdir;

foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="2">
my $maxThread = 10;
my $FCCutoffAry_ref = [0]; #----loop through each comparison and read the results
my $FPKMCutoff = 1; #---if both ref and qry are smaller than this value, it won't be tested in any enrichment test
my $refinePVal = 0.05;
my $FDRCutOffOut = 0.05;
my $maxGeneSetSize = 150;
my $minGeneSetSize = 10;
my $FUNCGOTermDBDir = '/Applications/Research/GeneOntology/FUNC/go_201209-termdb-tables/';
my $geneURLPrefix = "http://amoebadb.org/amoeba/showRecord.do?name=GeneRecordClasses.GeneRecordClass&source_id=<REPLACE_THIS>";
my $DEResultColumnIndexHsh_ref = {
	'AvgFpkm_qry'=>4, #---used for printing final list only, not functional
	'AvgFpkm_ref'=>5, #---used for printing final list only, not functional
	'DESeq1_linearFC'=>13, #---used for sorting FC cutoffs and subselect up/dn
	'DESeq1_modLog2FC'=>14, #---used for GSEA
	'DESeq1_qValue'=>16, #---used for printing final list only, not functional
	'significantDiffExp'=>27, #---used to define whether a gene is significantlt differentially expressed, any non-zero integer or string indicates a yes

#---For roman's data
#	
#	'AvgFpkm_qry'=>1, #---used for printing final list only, not functional
#	'AvgFpkm_ref'=>2, #---used for printing final list only, not functional
#	'DESeq1_linearFC'=>3, #---used for sorting FC cutoffs and subselect up/dn
#	'DESeq1_modLog2FC'=>4, #---used for GSEA
#	'DESeq1_qValue'=>5, #---used for printing final list only, not functional
#	'significantDiffExp'=>6, #---used to define whether a gene is significantlt differentially expressed, any non-zero integer or string indicates a yes

};
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineHardCodedFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedFilePath" num="3">
my $masterGAFPath = '/Volumes/A_MPro2TB/softwareForNGS/resources/genome/AmoebaDB/inUse/AmoebaDB-3.0/EHI/GOterms/gaf/AmoebaDB_v3.0.GOTerms.gaf';#---path of the GO association file of the genome, used to extract go terms of the gene list;
my $regularGOOboPath = '/Volumes/A_MPro2TB/softwareForNGS/resources/geneOntology/regular/gene_ontology.1_2.obo';#---path of the GO definition file of the whole gene ontology;
my $GOSlimOboPathAry_ref = [
	'/Volumes/A_MPro2TB/softwareForNGS/resources/geneOntology/map2Slim/pir.obo',
	'/Volumes/A_MPro2TB/softwareForNGS/resources/geneOntology/map2Slim/yeast.obo',
	'/Volumes/A_MPro2TB/softwareForNGS/resources/geneOntology/map2Slim/generic.obo'
]; #---path of the GO definition file of the GO Slim terms; can be more than one, one GO Slim per folder;

my $geneSetPathHsh_ref = {
	'GO'=>'/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GeneOntology/multiCuffDiffDeSeqRunnerResultAnanlyzer/v0.2/input/geneSets/GO.gmt',
	'KEGG'=>'/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GeneOntology/multiCuffDiffDeSeqRunnerResultAnanlyzer/v0.2/input/geneSets/KEGG.gmt',
	'FAM'=>'/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GeneOntology/multiCuffDiffDeSeqRunnerResultAnanlyzer/v0.2/input/geneSets/FAM.gmt',
	'INTERPRO'=>'/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/GeneOntology/multiCuffDiffDeSeqRunnerResultAnanlyzer/v0.2/input/geneSets/INTERPRO.gmt',
};

my $GSEABinPath= '/Applications/Research/GeneOntology/GSEA/commandLine/GSEA2-2.07/bin/gsea2-2.07.jar'; #---path of the GSEA jar, e.g. /Applications/Research/GeneOntology/GSEA/commandLine/GSEA2-2.07/bin/gsea2-2.07.jar
my $FUNCTermPath = '/Applications/Research/GeneOntology/FUNC/GeneToGOResources/term.txt';#---path of term.txt for FUNCGeneToGOScript
my $FUNCGraphPath = '/Applications/Research/GeneOntology/FUNC/GeneToGOResources/graph_path.txt';#---path for graph_path.txt for FUNCGeneToGOScript
my $FUNCGeneToGOScript = '/Applications/Research/GeneOntology/FUNC/func-0.4.5/src/tools/input_to_go_genes.pl';#---path of input_to_go_genes.pl to get the genes-to-go association
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: readGODefObo|931, readGeneDEResultInComparison|970, readGeneDEResultPathFile|1042, readGeneGAF|1073, readGeneInfo|1109
#	secondaryDependOnSub: reportStatus|1595
#
#<section ID="processInputData" num="4">
#---get all common info first
my ($rglrGODefOboHsh_ref, undef) = &readGODefObo($regularGOOboPath);#->931
my ($geneGAFInfoHsh_ref, $geneGAFOriginalLineHsh_ref, $geneGOProcessHsh_ref) = &readGeneGAF($masterGAFPath);#->1073
my ($geneInfoHsh_ref) = &readGeneInfo($geneInfoPath);#->1109
my ($DEComparisonInfoHsh_ref) = &readGeneDEResultPathFile($geneDEResultPathFile);#->1042
my ($sgnfcntDEGeneListHsh_ref, $DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref) = &readGeneDEResultInComparison($DEComparisonInfoHsh_ref, $FCCutoffAry_ref , $DEResultColumnIndexHsh_ref, $geneInfoHsh_ref);#->970
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_createEmptyVarToCollectDataFromSections
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="createEmptyVarToCollectDataFromSections" num="5">
my $geneBasedEnrichmentHsh_ref = {};
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_FUNCHyper
#	primaryDependOnSub: generateHyperFUNCMasterHtml|667, refineAndParseFUNCHyperResults|1178, runFUNCHyperInAllComparison|1852
#	secondaryDependOnSub: generateFUNCHyperInput|298, generateGeneListHtmlLinkForSVG|536, generateRamiGOTreesvg|750, generate_random_string|790, reportStatus|1595
#
#<section ID="FUNCHyper" num="6">
#---Run FUNC for all comparisons;
my ($FUNCHyperRunInfoHsh_ref) = &runFUNCHyperInAllComparison($sgnfcntDEGeneListHsh_ref, $FUNCHyperDir, $geneGAFInfoHsh_ref, $FUNCGOTermDBDir, $FUNCTermPath, $FUNCGraphPath, $FUNCGeneToGOScript, $FPKMCutoff, $DEResultDataHsh_ref, $maxThread);#->1852
my ($FUNCHyperResultPathHsh_ref) = &refineAndParseFUNCHyperResults($FUNCHyperRunInfoHsh_ref, $refinePVal, $FDRCutOffOut, $geneInfoHsh_ref, $DEResultDataHsh_ref, $sgnfcntDEGeneListHsh_ref, $FUNCHyperDir, $geneBasedEnrichmentHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEComparisonInfoHsh_ref);#->1178
my ($hyperFUNCMasterHtml) = &generateHyperFUNCMasterHtml($FUNCHyperResultPathHsh_ref, $dirForFUNCMasterHtml, $FDRCutOffOut, $refinePVal);#->667
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_GSEA
#	primaryDependOnSub: generateGSEAMasterHtml|486, refineAndParseGSEAResults|1427, runGSEAInAllComparison|1917
#	secondaryDependOnSub: generateGSEAInput|409, generateGeneListHtmlLinkForSVG|536, generateRamiGOTreesvg|750, reportStatus|1595
#
#<section ID="GSEA" num="7">
#---- run GSEA for all comparisons
my ($GSEARunInfoHsh_ref) = &runGSEAInAllComparison($DEResultDataHsh_ref, $GSEADir, $GSEABinPath, $maxGeneSetSize, $minGeneSetSize, $FPKMCutoff, $geneSetPathHsh_ref, $maxThread);#->1917
my ($GSEAResultDataHsh_ref, $GSEAResultPathHsh_ref, $GSEAGOSVGPathHsh_ref) = &refineAndParseGSEAResults($GSEARunInfoHsh_ref, $geneBasedEnrichmentHsh_ref, $DEResultDataHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEComparisonInfoHsh_ref);#->1427
my ($GSEAMasterHtmlPath) = &generateGSEAMasterHtml($GSEAResultPathHsh_ref, $GSEADir, $GSEAGOSVGPathHsh_ref);#->486
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 8_map2Slim
#	primaryDependOnSub: generateGOSlimMasterHtml|363, runMap2SlimInAllComparison|1970
#	secondaryDependOnSub: generateRamiGOTreesvg|750, printGeneListGAF|902, reportStatus|1595, runAndParseMap2Slim|1616
#
#<section ID="map2Slim" num="8">
#---- run map2Slim for all comparisons
my ($GOSlimResultPathHsh_ref, $GOSlimRootDir) = &runMap2SlimInAllComparison($sgnfcntDEGeneListHsh_ref, $DEResultDataHsh_ref, $geneGAFOriginalLineHsh_ref, $GOSlimOboPathAry_ref, $regularGOOboPath, $fullDEInfoHeaderAry_ref, $outDir);#->1970
my $dirForGOSlimMasterHtml = $GOSlimRootDir;
my ($GOSlimMasterHtml) = &generateGOSlimMasterHtml($GOSlimResultPathHsh_ref, $dirForGOSlimMasterHtml);#->363
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 9_outputData
#	primaryDependOnSub: generateGeneSummaryMasterHtml|632, generateMasterIndexHtml|710, printGeneBasedSummary|851
#	secondaryDependOnSub: reportStatus|1595
#
#<section ID="outputData" num="9">
#---- print gene based summary
my $geneBasedSummaryPathHsh_ref =  &printGeneBasedSummary($DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref, $geneBasedEnrichmentHsh_ref, $geneBasedSummaryOutdir);#->851
my $dirForGeneSummaryMasterHtml = $geneBasedSummaryOutdir;
my $geneSummaryMasterHtmlPath = &generateGeneSummaryMasterHtml($geneBasedSummaryPathHsh_ref, $dirForGeneSummaryMasterHtml);#->632
#---- generation of mater index
&generateMasterIndexHtml($hyperFUNCMasterHtml, $GSEAMasterHtmlPath, $GOSlimMasterHtml, $geneSummaryMasterHtmlPath, $outDir);#->710
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 10_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|818
#	secondaryDependOnSub: currentTime|280
#
#<section ID="finishingTasks" num="10">
&printCMDLogOrFinishMessage("finishMessage");#->818
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	GeneOntology [n=9]:
#		generateFUNCHyperInput, generateGSEAInput, generateRamiGOTreesvg
#		printGeneListGAF, readGODefObo, readGeneGAF
#		refineAndParseFUNCHyperResults, refineAndParseGSEAResults, runAndParseMap2Slim
#
#	HTML [n=6]:
#		generateGOSlimMasterHtml, generateGSEAMasterHtml, generateGeneListHtmlLinkForSVG
#		generateGeneSummaryMasterHtml, generateHyperFUNCMasterHtml, generateMasterIndexHtml
#
#	general [n=6]:
#		currentTime, generate_random_string, printCMDLogOrFinishMessage
#		readParameters, reportStatus, runFUNCHyperInAllComparison
#
#	getTextInfo [n=1]:
#		readGeneInfo
#
#	specific [n=5]:
#		printGeneBasedSummary, readGeneDEResultInComparison, readGeneDEResultPathFile
#		runGSEAInAllComparison, runMap2SlimInAllComparison
#
#====================================================================================================================================================#

sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|818, reportStatus|1595
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|66, 10_finishingTasks|240
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 55, 838, 841, 846, 1611
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateFUNCHyperInput {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: >none
#	appearInSub: runFUNCHyperInAllComparison|1852
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_FUNCHyper|186
#	input: $DEResultDataHsh_ref, $FPKMCutoff, $FUNCGOTermDBDir, $FUNCGeneToGOScript, $FUNCGraphPath, $FUNCHyperDir, $FUNCTermPath, $allUpDn, $comparison, $cutoff, $geneGAFInfoHsh_ref, $sgnfcntDEGeneListHsh_ref
#	output: $FUNCHyperOutDir, $funcGOToGeneCMD, $funcGOToGeneOutputPath, $runFUNCHyperCmd
#	toCall: my ($runFUNCHyperCmd, $FUNCHyperOutDir, $funcGOToGeneOutputPath, $funcGOToGeneCMD) = &generateFUNCHyperInput($FUNCHyperDir, $sgnfcntDEGeneListHsh_ref, $geneGAFInfoHsh_ref, $comparison, $cutoff, $FUNCGOTermDBDir, $FUNCTermPath, $FUNCGraphPath, $FUNCGeneToGOScript, $allUpDn, $FPKMCutoff, $DEResultDataHsh_ref);
#	calledInLine: 1884
#....................................................................................................................................................#

	my ($FUNCHyperDir, $sgnfcntDEGeneListHsh_ref, $geneGAFInfoHsh_ref, $comparison, $cutoff, $FUNCGOTermDBDir, $FUNCTermPath, $FUNCGraphPath, $FUNCGeneToGOScript, $allUpDn, $FPKMCutoff, $DEResultDataHsh_ref) = @_;

	my $FUNCHyperOutDir = "$FUNCHyperDir/$comparison/$cutoff/$allUpDn/output";
	my $FUNCHyperInDir = "$FUNCHyperDir/$comparison/$cutoff/$allUpDn/input";

	system ("mkdir -p -m 777 $FUNCHyperOutDir");
	system ("mkdir -p -m 777 $FUNCHyperInDir");
	
	#---generate the chip from DEResultDataHsh
	my $FUNCHyperTSVPath = "$FUNCHyperInDir/FUNCHyperTSV.$allUpDn.tsv";
	open (TSV, ">$FUNCHyperTSVPath");
	foreach my $geneID (sort {$a cmp $b} keys %{$geneGAFInfoHsh_ref}) {
		
		#----process genes only inside the scope of DE
		if ($DEResultDataHsh_ref->{$comparison}{$geneID}) {
			my $AvgFpkm_qry = $DEResultDataHsh_ref->{$comparison}{$geneID}{'AvgFpkm_qry'};
			my $AvgFpkm_ref = $DEResultDataHsh_ref->{$comparison}{$geneID}{'AvgFpkm_ref'};
			next if $AvgFpkm_qry < $FPKMCutoff and $AvgFpkm_ref < $FPKMCutoff;
			
			my $significantDiffExp = 0;
			if (exists $sgnfcntDEGeneListHsh_ref->{$comparison}{$cutoff}{$geneID}) {
				my $DESeq1_linearFC = $sgnfcntDEGeneListHsh_ref->{$comparison}{$cutoff}{$geneID};
				$significantDiffExp = 1 
				if (($DESeq1_linearFC > 0 and $allUpDn eq 'up')
				or	($DESeq1_linearFC < 0 and $allUpDn eq 'dn')
				or	($allUpDn eq 'all'));
			}
			foreach my $GOID (sort {$a cmp $b} keys %{$geneGAFInfoHsh_ref->{$geneID}}) {
				print TSV join '', ((join "\t", ($geneID, $GOID, $significantDiffExp)), "\n");
			}
		}
	}
	close TSV;
	
	my $funcGOToGeneOutputPath = "$FUNCHyperOutDir/funcGOToGene.log.tsv";
	my $funcGOToGeneCMD = "perl $FUNCGeneToGOScript $FUNCTermPath $FUNCGraphPath $FUNCHyperTSVPath >$funcGOToGeneOutputPath";
	
	my $inputTSV = "-i $FUNCHyperTSVPath";
	my $outputDir = "-o $FUNCHyperOutDir";
	my $GOTable = "-t $FUNCGOTermDBDir";
	my $FUNCCutoff = "-c 1";
	my $randSet = "-r 1000";
	
	my $runFUNCHyperCmd = join " ", ('func_hyper', $inputTSV, $GOTable, $outputDir, $FUNCCutoff, $randSet);
	my $cmdPath = "$FUNCHyperInDir/runFUNCHyperCmd.txt";
	open (CMD, ">$cmdPath");
	print CMD $runFUNCHyperCmd."\n";
	print CMD $funcGOToGeneCMD."\n";
	close CMD;
	
	return $runFUNCHyperCmd, $FUNCHyperOutDir, $funcGOToGeneOutputPath, $funcGOToGeneCMD;
}
sub generateGOSlimMasterHtml {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 8_map2Slim|212
#	secondaryAppearInSection: >none
#	input: $GOSlimResultPathHsh_ref, $dirForGOSlimMasterHtml
#	output: $GOSlimMasterHtml
#	toCall: my ($GOSlimMasterHtml) = &generateGOSlimMasterHtml($GOSlimResultPathHsh_ref, $dirForGOSlimMasterHtml);
#	calledInLine: 220
#....................................................................................................................................................#

	my ($GOSlimResultPathHsh_ref, $dirForGOSlimMasterHtml) = @_;
	my $GOSlimMasterHtml = "$dirForGOSlimMasterHtml/master_GOSlim.html";
	
	open (MASTERGOSLIMHTML, ">", $GOSlimMasterHtml);
	print MASTERGOSLIMHTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print MASTERGOSLIMHTML join "", ('<html>', "\n");
	print MASTERGOSLIMHTML join "", ("<head><title>GO Slim summary of differentially expressed genes</title></head>\n");
	print MASTERGOSLIMHTML join "", ('<body>', "\n");
	print MASTERGOSLIMHTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");

	foreach my $cutoff (sort keys %{$GOSlimResultPathHsh_ref}) {
		print MASTERGOSLIMHTML join "", ('</div><div>', "\n");
		print MASTERGOSLIMHTML join "", ("<h4>GO Slim summary of differentially expressed genes at cutoff $cutoff</h4>\n");
		print MASTERGOSLIMHTML join "", ('<ul>', "\n");
		foreach my $comparison (sort keys %{$GOSlimResultPathHsh_ref->{$cutoff}}) {
			
			my $svgPath = $GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'svgPath'};
			my $geneBasedOneLinerGOSlimResultPath = $GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'geneBasedOneLinerGOSlimResultPath'};
			my $GOBasedOneLinerGOSlimResultPath = $GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'GOBasedOneLinerGOSlimResultPath'};
			my $geneGOSlimOneToOneResultPath = $GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'geneGOSlimOneToOneResultPath'};

			$svgPath =~ s/$dirForGOSlimMasterHtml//;$svgPath =~ s/^\///;
			$geneBasedOneLinerGOSlimResultPath =~ s/$dirForGOSlimMasterHtml//;$geneBasedOneLinerGOSlimResultPath =~ s/^\///;
			$GOBasedOneLinerGOSlimResultPath =~ s/$dirForGOSlimMasterHtml//;$GOBasedOneLinerGOSlimResultPath =~ s/^\///;
			$geneGOSlimOneToOneResultPath =~ s/$dirForGOSlimMasterHtml//;$geneGOSlimOneToOneResultPath =~ s/^\///;

			print MASTERGOSLIMHTML "<li>$comparison<a href=\'$svgPath\'>&nbsp[GOTree]&nbsp</a><a href=\'$GOBasedOneLinerGOSlimResultPath\'>&nbsp[GOBasedTable]&nbsp</a><a href=\'$geneBasedOneLinerGOSlimResultPath\'>&nbsp[GeneBasedTable]&nbsp</a><a href=\'$geneGOSlimOneToOneResultPath\'>&nbsp[OneGeneToOneGOTable]&nbsp</a>\n";
		}
	}
	close MASTERGOSLIMHTML;
	
	return $GOSlimMasterHtml;
}
sub generateGSEAInput {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: >none
#	appearInSub: runGSEAInAllComparison|1917
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 7_GSEA|199
#	input: $DEResultDataHsh_ref, $FPKMCutoff, $GSEABinPath, $GSEADir, $comparison, $geneSetPathHsh_ref, $maxGeneSetSize, $minGeneSetSize
#	output: $GSEACmdAndDirHsh_ref
#	toCall: my ($GSEACmdAndDirHsh_ref) = &generateGSEAInput($GSEADir, $DEResultDataHsh_ref, $geneSetPathHsh_ref, $comparison, $GSEABinPath, $maxGeneSetSize, $minGeneSetSize, $FPKMCutoff);
#	calledInLine: 1941
#....................................................................................................................................................#

	my ($GSEADir, $DEResultDataHsh_ref, $geneSetPathHsh_ref, $comparison, $GSEABinPath, $maxGeneSetSize, $minGeneSetSize, $FPKMCutoff) = @_;

	system ("mkdir -p -m 777 $GSEADir/$comparison/input/");
	system ("mkdir -p -m 777 $GSEADir/$comparison/output/");
	
	#---generate the chip from DEResultDataHsh
	my $chipPath = "$GSEADir/$comparison/input/geneAnnotation.chip";
	my $rnkPath = "$GSEADir/$comparison/input/DESeq1_modLog2FC.rnk";
	open (CHIP, ">$chipPath");
	open (RNK, ">$rnkPath");
	print CHIP join '', ((join "\t", ('Probe Set ID', 'Gene Symbol', 'Gene Title')), "\n");
	foreach my $geneID (sort {$a cmp $b} keys %{$DEResultDataHsh_ref->{$comparison}}) {
		my $AvgFpkm_ref = $DEResultDataHsh_ref->{$comparison}{$geneID}{'AvgFpkm_ref'};
		my $AvgFpkm_qry = $DEResultDataHsh_ref->{$comparison}{$geneID}{'AvgFpkm_qry'};  
		next if not $AvgFpkm_ref or not $AvgFpkm_qry; #---genes have annotation but do not exists in the dataset
		next if $AvgFpkm_ref < $FPKMCutoff and $AvgFpkm_qry < $FPKMCutoff;

		my $DESeq1_modLog2FC = $DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_modLog2FC'};
		my $geneName = $DEResultDataHsh_ref->{$comparison}{$geneID}{'geneName'};
		print CHIP join '', ((join "\t", ($geneID, $geneID, $geneName)), "\n");
		print RNK join '', ((join "\t", ($geneID, $DESeq1_modLog2FC)), "\n");

	}
	close CHIP;
	close RNK;
	
	my $GSEACmdAndDirHsh_ref = {};
	
	#----define the GSEA parameters
	my $cmdPath = "$GSEADir/$comparison/input/runGSEACmd.txt";
	open (CMD, ">$cmdPath");
	foreach my $geneSetType (keys %{$geneSetPathHsh_ref}) {
		my $geneSetPath = $geneSetPathHsh_ref->{$geneSetType};
		my $GSEAOutDir = "$GSEADir/$comparison/output/$geneSetType/";
		system ("mkdir -p -m 777 $GSEAOutDir");
		my $cp = "-cp $GSEABinPath";
		my $Xmx512m = "-Xmx512m xtools.gsea.GseaPreranked";
		my $chip = "-chip $chipPath";
		my $collapse = "-collapse false";
		my $mode = "-mode Max_probe";
		my $norm = "-norm meandiv";
		my $nperm = "-nperm 1000";
		my $rnk = "-rnk $rnkPath";
		my $scoring_scheme = "-scoring_scheme weighted";
		my $rpt_label = "-rpt_label $comparison";
		my $include_only_symbols = "-include_only_symbols true";
		my $make_sets = "-make_sets true";
		my $plot_top_x = "-plot_top_x 500";
		my $rnd_seed = "-rnd_seed ";
		my $set_max = "-set_max $maxGeneSetSize";
		my $set_min = "-set_min $minGeneSetSize";
		my $zip_report = "-zip_report false";
		my $gui = "-gui false";
		my $gmx = "-gmx $geneSetPath";
		my $out = "-out $GSEAOutDir";
		my $runGSEACmd = join " ", ('java', $cp, $Xmx512m, $gmx, $chip, $collapse, $mode, $norm, $nperm, $rnk, $scoring_scheme, $rpt_label, $include_only_symbols, $make_sets, $plot_top_x, $rnd_seed, $set_max, $set_min, $zip_report, $gui, $out);
		$GSEACmdAndDirHsh_ref->{$geneSetType}{'runGSEACmd'} = $runGSEACmd;
		$GSEACmdAndDirHsh_ref->{$geneSetType}{'GSEAOutDir'} = $GSEAOutDir;
		print CMD $runGSEACmd."\n";
	}
	close CMD;
	
	return ($GSEACmdAndDirHsh_ref);
}
sub generateGSEAMasterHtml {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 7_GSEA|199
#	secondaryAppearInSection: >none
#	input: $GSEAGOSVGPathHsh_ref, $GSEAResultPathHsh_ref, $dirForGSEAMasterHtml
#	output: $GSEAMasterHtmlPath
#	toCall: my ($GSEAMasterHtmlPath) = &generateGSEAMasterHtml($GSEAResultPathHsh_ref, $dirForGSEAMasterHtml, $GSEAGOSVGPathHsh_ref);
#	calledInLine: 207
#....................................................................................................................................................#

	my ($GSEAResultPathHsh_ref, $dirForGSEAMasterHtml, $GSEAGOSVGPathHsh_ref) = @_;
	my $GSEAMasterHtmlPath = "$dirForGSEAMasterHtml/master_GSEA.html";
	
	open (MASTERGSEAHTML, ">", $GSEAMasterHtmlPath);
	print MASTERGSEAHTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print MASTERGSEAHTML join "", ('<html>', "\n");
	print MASTERGSEAHTML join "", ('<head><title>Master index for all GSEA</title></head>',"\n");
	print MASTERGSEAHTML join "", ('<body>', "\n");
	print MASTERGSEAHTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	foreach my $geneSetType (sort keys %{$GSEAResultPathHsh_ref}) {
		print MASTERGSEAHTML join "", ('</div><div>', "\n");
		print MASTERGSEAHTML join "", ('<h4>GSEA Enrichment of ', $geneSetType, ' genesets</h4>', "\n");
		print MASTERGSEAHTML join "", ('<ul>', "\n");
		foreach my $comparison (sort keys %{$GSEAResultPathHsh_ref->{$geneSetType}}) {
			my $indexHtmlRelative = $GSEAResultPathHsh_ref->{$geneSetType}{$comparison}{'idx'};
			$indexHtmlRelative =~ s/$dirForGSEAMasterHtml//;
			$indexHtmlRelative =~ s/^\///;
			print MASTERGSEAHTML "<li><a href=\'$indexHtmlRelative\'>$comparison [mainIndex]</a>";
			if ($geneSetType eq "GO") {
				foreach my $FDRCutOff (sort keys %{$GSEAGOSVGPathHsh_ref->{$comparison}}) {
					my $GSEAGOSVGRelative = $GSEAGOSVGPathHsh_ref->{$comparison}{$FDRCutOff};
					$GSEAGOSVGRelative =~ s/$dirForGSEAMasterHtml//;
					$GSEAGOSVGRelative =~ s/^\///;
					print MASTERGSEAHTML "&nbsp&nbsp&nbsp<a href=\'$GSEAGOSVGRelative\'>GoTree[FDR$FDRCutOff]</a>";
				}
			}
			print MASTERGSEAHTML "\n";
		}
		print MASTERGSEAHTML join "", ('</ul>', "\n");
	}
	
	print MASTERGSEAHTML join "", ('</html>', "\n");
	print MASTERGSEAHTML join "", ('</body>', "\n");
	close MASTERGSEAHTML;
	
	return $GSEAMasterHtmlPath;
}
sub generateGeneListHtmlLinkForSVG {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: refineAndParseFUNCHyperResults|1178, refineAndParseGSEAResults|1427
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_FUNCHyper|186, 7_GSEA|199
#	input: $DEResultDataHsh_ref, $analysisType, $comparison, $geneIDInGOGeneSetHsh_ref, $geneURLPrefix, $qrySample, $refSample, $rglrGODefOboHsh_ref, $svgDir, $svgName, $svgPath
#	output: none
#	toCall: &generateGeneListHtmlLinkForSVG($geneIDInGOGeneSetHsh_ref, $svgPath, $svgName, $svgDir, $comparison, $refSample, $qrySample, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEResultDataHsh_ref, $analysisType);
#	calledInLine: 1411, 1584
#....................................................................................................................................................#

	my ($geneIDInGOGeneSetHsh_ref, $svgPath, $svgName, $svgDir, $comparison, $refSample, $qrySample, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEResultDataHsh_ref, $analysisType) = @_;
	
	open (svgIN, "<", $svgPath);
	my @svgBufferAry = <svgIN>;
	close svgIN;
	my @svgOutAry;
	
	my $dirForHtmlLinks = $svgDir."/GOGeneHtml/$svgName/";
	system ("mkdir -p -m 777 $dirForHtmlLinks");
	for (my $i=0; $i<=$#svgBufferAry; $i++) {
		if ($svgBufferAry[$i] =~ m/>GO:(\d+)</) {
			my $GOID = "GO:".$1;
			my $GOLINK = "GO_".$1;
			if (exists $geneIDInGOGeneSetHsh_ref->{$GOID}) {
				my $htmlLink = "GOGeneHtml/$svgName/$GOLINK.html";
				my $TSVLink = "GOGeneHtml/$svgName/$GOLINK.tsv";
				$svgBufferAry[$i] = "<a xlink:href=\'".$htmlLink."\'>".$svgBufferAry[$i];
				$svgBufferAry[$i+1] = $svgBufferAry[$i+1]."</a>";
				$i++;

				my $GOName = my $GODef = my $upOrDn = "obsolete"; 
				$GOName =  $rglrGODefOboHsh_ref->{$GOID}{"name"} if exists $rglrGODefOboHsh_ref->{$GOID}{"name"};
				$GODef =  $rglrGODefOboHsh_ref->{$GOID}{"def"} if exists $rglrGODefOboHsh_ref->{$GOID}{"def"};
				$upOrDn = $geneIDInGOGeneSetHsh_ref->{$GOID}{'upOrDn'} if exists $geneIDInGOGeneSetHsh_ref->{$GOID}{'upOrDn'};
				$GODef =~ s/\[.+\]//g;
				
				open (HTMLLINK, ">", "$svgDir/$htmlLink");
				open (TSVLINK, ">", "$svgDir/$TSVLink");
				print HTMLLINK "<html><head><title>$GOID [$GOName] $upOrDn regulated in $qrySample of $comparison</title></head><body>\n";
				print HTMLLINK "<div class='richTable'>\n";
				print HTMLLINK "<table cols='11' border='1'>\n";
				print HTMLLINK "<caption class='table'>";
				print HTMLLINK "<b>$upOrDn</b> regulated in $qrySample of $comparison in $analysisType analysis<br>";
				print HTMLLINK "<a href=\'http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=$GOID\'>$GOID</a> <b>[$GOName]</b> $GODef<br>";
				print HTMLLINK "<a href=\"$GOLINK.tsv\">click here to open in excel</a>";
				print HTMLLINK "</caption>\n";
				print HTMLLINK "<th class='richTable'></th>\n";
				print HTMLLINK "<th class='richTable'>GeneID</th>\n";
				print HTMLLINK "<th class='richTable'>GeneName</th>\n";
				print HTMLLINK "<th class='richTable'>DESeq1_linearFC</th>\n";
				print HTMLLINK "<th class='richTable'>fpkm in $refSample</th>\n";
				print HTMLLINK "<th class='richTable'>fpkm in $qrySample</th>\n";
				print HTMLLINK "<th class='richTable'>core(in GSEA) or significantly differential(in others)</th>\n";
				
				print TSVLINK "$GOID [$GOName] $upOrDn regulated in $qrySample of $comparison\n";
				print TSVLINK "$GODef\n";
				print TSVLINK join "", ((join "\t", ('GeneID', 'GeneName', 'DESeq1_linearFC', "fpkm in $refSample", "fpkm in $qrySample", 'core(in GSEA) or significantly differential(in others)')), "\n");
				
				
				my $geneNum = 0;
				foreach my $geneID (sort {$DEResultDataHsh_ref->{$comparison}{$b}{'DESeq1_linearFC'} <=> $DEResultDataHsh_ref->{$comparison}{$a}{'DESeq1_linearFC'}} keys %{$geneIDInGOGeneSetHsh_ref->{$GOID}{'gene'}}) {
					$geneNum++;
					my $AvgFpkm_ref = $DEResultDataHsh_ref->{$comparison}{$geneID}{'AvgFpkm_ref'};
					my $AvgFpkm_qry = $DEResultDataHsh_ref->{$comparison}{$geneID}{'AvgFpkm_qry'};  
					my $DESeq1_linearFC = $DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_linearFC'};
					my $geneName = $DEResultDataHsh_ref->{$comparison}{$geneID}{'geneName'};
					my $geneURL = $geneURLPrefix;
					my $coreEnrichment = $geneIDInGOGeneSetHsh_ref->{$GOID}{'gene'}{$geneID};
					$geneURL =~ s/<REPLACE_THIS>/$geneID/;
					print HTMLLINK "<tr>\n";
					print HTMLLINK "<td class='lessen'>$geneNum</td>\n";
					print HTMLLINK "<td><a href=\'$geneURL\'>$geneID</a></td>\n";
					print HTMLLINK "<td>$geneName</td>\n";
					print HTMLLINK "<td>$DESeq1_linearFC</td>\n";
					print HTMLLINK "<td>$AvgFpkm_ref</td>\n";
					print HTMLLINK "<td>$AvgFpkm_qry</td>\n";
					print HTMLLINK "<td>$coreEnrichment</td>\n";
					print TSVLINK join "", ((join "\t", ($geneID, $geneName, $DESeq1_linearFC, $AvgFpkm_ref, $AvgFpkm_qry, $coreEnrichment)), "\n");
				}
				print HTMLLINK "</table></div></body></html>\n";
				close HTMLLINK;
				close TSVLINK;
			}
		}
	}
	
	open (svgOUT, ">", $svgPath);
	foreach (@svgBufferAry) {
		print svgOUT $_;
	}
	close svgOUT;
	
}
sub generateGeneSummaryMasterHtml {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 9_outputData|225
#	secondaryAppearInSection: >none
#	input: $dirForGeneSummaryMasterHtml, $geneBasedSummaryPathHsh_ref
#	output: $geneSummaryMasterHtmlPath
#	toCall: my ($geneSummaryMasterHtmlPath) = &generateGeneSummaryMasterHtml($geneBasedSummaryPathHsh_ref, $dirForGeneSummaryMasterHtml);
#	calledInLine: 233
#....................................................................................................................................................#

	my ($geneBasedSummaryPathHsh_ref, $dirForGeneSummaryMasterHtml) = @_;
	my $geneSummaryMasterHtmlPath = "$dirForGeneSummaryMasterHtml/master_geneSummary.html";

	open (MASTERGENESUMHTML, ">", $geneSummaryMasterHtmlPath);
	print MASTERGENESUMHTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print MASTERGENESUMHTML join "", ('<html>', "\n");
	print MASTERGENESUMHTML join "", ("<head><title>Gene based summary of all analyses</title></head>\n");
	print MASTERGENESUMHTML join "", ('<body>', "\n");
	print MASTERGENESUMHTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");
	print MASTERGENESUMHTML join "", ('</div><div>', "\n");
	print MASTERGENESUMHTML join "", ("<h4>Gene based summary of all analyses</h4>\n");
	print MASTERGENESUMHTML join "", ('<ul>', "\n");
	foreach my $comparison (sort keys %{$geneBasedSummaryPathHsh_ref}) {
			
		my $geneBasedSummaryPath = $geneBasedSummaryPathHsh_ref->{$comparison};
		$geneBasedSummaryPath =~ s/$dirForGeneSummaryMasterHtml//;$geneBasedSummaryPath =~ s/^\///;
		print MASTERGENESUMHTML "<li><a href=\'$geneBasedSummaryPath\'>$comparison</a>\n";
	}
	close MASTERGENESUMHTML;

	return $geneSummaryMasterHtmlPath;
}
sub generateHyperFUNCMasterHtml {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 6_FUNCHyper|186
#	secondaryAppearInSection: >none
#	input: $FDRCutOffOut, $FUNCHyperResultPathHsh_ref, $dirForHyperFUNCMasterHtml, $refinePVal
#	output: $hyperFUNCMasterHtml
#	toCall: my ($hyperFUNCMasterHtml) = &generateHyperFUNCMasterHtml($FUNCHyperResultPathHsh_ref, $dirForHyperFUNCMasterHtml, $FDRCutOffOut, $refinePVal);
#	calledInLine: 194
#....................................................................................................................................................#

	my ($FUNCHyperResultPathHsh_ref, $dirForHyperFUNCMasterHtml, $FDRCutOffOut, $refinePVal) = @_;
	my $hyperFUNCMasterHtml = "$dirForHyperFUNCMasterHtml/master_hyperFUNC.html";
	
	open (MASTERFUNCHTML, ">", $hyperFUNCMasterHtml);
	print MASTERFUNCHTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print MASTERFUNCHTML join "", ('<html>', "\n");
	print MASTERFUNCHTML join "", ("<head><title>Master index for all hyperFUNC GO enrichment analyses at FDR $FDRCutOffOut and refinement pval $refinePVal</title></head>\n");
	print MASTERFUNCHTML join "", ('<body>', "\n");
	print MASTERFUNCHTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");

	foreach my $cutoff (sort keys %{$FUNCHyperResultPathHsh_ref}) {
		print MASTERFUNCHTML join "", ('</div><div>', "\n");
		print MASTERFUNCHTML join "", ("<h4>hyperFUNC GO Enrichment of differential expressed gene using cutoff at $cutoff</h4>\n");
		print MASTERFUNCHTML join "", ('<ul>', "\n");
		foreach my $comparison (sort keys %{$FUNCHyperResultPathHsh_ref->{$cutoff}}) {
			
			my $GOTsvPath = $FUNCHyperResultPathHsh_ref->{$cutoff}{$comparison}{'GOTsvPath'};
			my $geneTsvPath = $FUNCHyperResultPathHsh_ref->{$cutoff}{$comparison}{'geneTsvPath'};
			my $svgPath = $FUNCHyperResultPathHsh_ref->{$cutoff}{$comparison}{'svgPath'};
			$GOTsvPath =~ s/$dirForHyperFUNCMasterHtml//;$GOTsvPath =~ s/^\///;
			$geneTsvPath =~ s/$dirForHyperFUNCMasterHtml//;$geneTsvPath =~ s/^\///;
			$svgPath =~ s/$dirForHyperFUNCMasterHtml//;$svgPath =~ s/^\///;
			print MASTERFUNCHTML "<li>$comparison<a href=\'$svgPath\'>&nbsp[GOTree]&nbsp</a><a href=\'$GOTsvPath\'>&nbsp[GOTable]&nbsp</a><a href=\'$geneTsvPath\'>&nbsp[GeneTable]&nbsp</a>\n";
		}
	}
	close MASTERFUNCHTML;
	
	return $hyperFUNCMasterHtml;
	
}
sub generateMasterIndexHtml {
#....................................................................................................................................................#
#	subroutineCategory: HTML
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 9_outputData|225
#	secondaryAppearInSection: >none
#	input: $GOSlimMasterHtml, $GSEAMasterHtmlPath, $geneSummaryMasterHtmlPath, $hyperFUNCMasterHtml, $outDir
#	output: none
#	toCall: &generateMasterIndexHtml($hyperFUNCMasterHtml, $GSEAMasterHtmlPath, $GOSlimMasterHtml, $geneSummaryMasterHtmlPath, $outDir);
#	calledInLine: 235
#....................................................................................................................................................#

	my ($hyperFUNCMasterHtml, $GSEAMasterHtmlPath, $GOSlimMasterHtml, $geneSummaryMasterHtmlPath, $outDir) = @_;
	
	my $indexMasterHtml = "$outDir/master_index.html";
	
	open (MASTERINDEXHTML, ">", $indexMasterHtml);
	print MASTERINDEXHTML join "", ('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/1999/REC-html401-19991224/loose.dtd">', "\n");
	print MASTERINDEXHTML join "", ('<html>', "\n");
	print MASTERINDEXHTML join "", ("<head><title>GO Slim summary of differentially expressed genes</title></head>\n");
	print MASTERINDEXHTML join "", ('<body>', "\n");
	print MASTERINDEXHTML join "", ('<div id="footer" style="width: 905; height: 35">', "\n");

	print MASTERINDEXHTML join "", ('</div><div>', "\n");
	print MASTERINDEXHTML join "", ("<h4>Click the following links to access the relevant information</h4>\n");
	print MASTERINDEXHTML join "", ('<ul>', "\n");

	$GOSlimMasterHtml =~ s/$outDir//;$GOSlimMasterHtml =~ s/^\///;
	$GSEAMasterHtmlPath =~ s/$outDir//;$GSEAMasterHtmlPath =~ s/^\///;
	$hyperFUNCMasterHtml =~ s/$outDir//;$hyperFUNCMasterHtml =~ s/^\///;
	$geneSummaryMasterHtmlPath =~ s/$outDir//;$geneSummaryMasterHtmlPath =~ s/^\///;

	print MASTERINDEXHTML "<li><a href=\'$GOSlimMasterHtml\'>Map2Slim Gene Ontology Summary of differentially expressed genes</a>\n";
	print MASTERINDEXHTML "<li><a href=\'$GSEAMasterHtmlPath\'>Gene Set Enrichment Analysis all gene expression profiles</a>\n";
	print MASTERINDEXHTML "<li><a href=\'$hyperFUNCMasterHtml\'>hyperFUNC analysis of enriched GO terms in differentially expressed genes</a>\n";
	print MASTERINDEXHTML "<li><a href=\'$geneSummaryMasterHtmlPath\'>Table summary expression data and functional enrichment data</a>\n";
	close MASTERINDEXHTML;

}
sub generateRamiGOTreesvg {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: reportStatus|1595
#	appearInSub: refineAndParseFUNCHyperResults|1178, refineAndParseGSEAResults|1427, runMap2SlimInAllComparison|1970
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_FUNCHyper|186, 7_GSEA|199, 8_map2Slim|212
#	input: $GOColorRamiGOHsh_ref, $RamiGOOutDir, $svgName
#	output: none
#	toCall: &generateRamiGOTreesvg($GOColorRamiGOHsh_ref, $RamiGOOutDir, $svgName);
#	calledInLine: 1407, 1580, 2012
#....................................................................................................................................................#

	my ($GOColorRamiGOHsh_ref, $RamiGOOutDir, $svgName) = @_;

	my $GOIDAry_ref = [];
	my $colorAry_ref = [];
	
	&reportStatus("Running RamiGO for $svgName", 10, "\n");#->1595

	foreach my $GOID (sort keys %{$GOColorRamiGOHsh_ref}) {
		push @{$GOIDAry_ref}, "\"$GOID\"";
		push @{$colorAry_ref}, "\"$GOColorRamiGOHsh_ref->{$GOID}\"";
	}
	my $GOIDStr = join ",", @{$GOIDAry_ref};
	my $colorStr = join ",", @{$colorAry_ref};
	
	system ("mkdir -p -m 777 $RamiGOOutDir/tmp  >/dev/null 2>&1");
	my $RamiGORScriptPath = "$RamiGOOutDir/tmp/$svgName.R";
	my $RamiGOsvgPath = "$RamiGOOutDir/$svgName.svg";
	open (RAMIGOR, ">", "$RamiGORScriptPath");
	print RAMIGOR "library(RamiGO);\n";
	print RAMIGOR "goIDs<-c($GOIDStr);\n";
	print RAMIGOR "color<-c($colorStr);\n";
	print RAMIGOR "svgRes <- getAmigoTree(goIDs=goIDs, color=color, filename=\"$RamiGOsvgPath\", picType=\"svg\", saveResult=TRUE);\n";
	close RAMIGOR;
	
	system ("R --slave --no-restore --vanilla --file=$RamiGORScriptPath >$RamiGOOutDir/tmp/R.error.log.txt 2>&1");
	
}
sub generate_random_string {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: refineAndParseFUNCHyperResults|1178
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 6_FUNCHyper|186
#	input: $length_of_randomstring
#	output: $random_string
#	toCall: my ($random_string) = &generate_random_string($length_of_randomstring);
#	calledInLine: 1302
#....................................................................................................................................................#

	#----code from http://th.atguy.com/mycode/&generate_random_string/
	my ($length_of_randomstring) = @_;# the length of 
	 # the random string to generate

	my @chars=('a'..'z','A'..'Z','0'..'9','_');
	my $random_string;
	foreach (1..$length_of_randomstring) 
	{
		# rand @chars will generate a random 
		# number between 0 and scalar @chars
		$random_string.=$chars[rand @chars];
	}

	return ($random_string);
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|280
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|66, 10_finishingTasks|240
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 72, 245
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->280
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->280
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->280
		print "=========================================================================\n\n";
	}
}
sub printGeneBasedSummary {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 9_outputData|225
#	secondaryAppearInSection: >none
#	input: $DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref, $geneBasedEnrichmentHsh_ref, $geneBasedSummaryOutdir
#	output: $geneBasedSummaryPathHsh_ref
#	toCall: my ($geneBasedSummaryPathHsh_ref) = &printGeneBasedSummary($DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref, $geneBasedEnrichmentHsh_ref, $geneBasedSummaryOutdir);
#	calledInLine: 231
#....................................................................................................................................................#

	my ($DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref, $geneBasedEnrichmentHsh_ref, $geneBasedSummaryOutdir) = @_;

	my $geneBasedSummaryPathHsh_ref = {};

	system ("mkdir -p -m 777 $geneBasedSummaryOutdir >/dev/null 2>&1");
	foreach my $comparison (sort keys %{$DEResultDataHsh_ref}) {

		my $geneBasedSummaryPath = "$geneBasedSummaryOutdir/$comparison.geneBasedSummary.tsv";
		$geneBasedSummaryPathHsh_ref->{$comparison} = $geneBasedSummaryPath;
		
		open GENEBASEDSUM, ">", "$geneBasedSummaryPath";
		my @outputAry = ('geneID', 'geneName');
		push @outputAry, @{$fullDEInfoHeaderAry_ref};
		foreach my $enrichmentInfoType (sort keys %{$geneBasedEnrichmentHsh_ref->{$comparison}}) {
			push @outputAry, $enrichmentInfoType;
		}
		print GENEBASEDSUM join "", (join "\t", @outputAry), "\n";
		
		foreach my $geneID (sort keys %{$DEResultDataHsh_ref->{$comparison}}) {
			if (not $DEResultDataHsh_ref->{$comparison}{$geneID}{'fullDEInfo'}) {
				&reportStatus("WARNING: $geneID is missing in DEResultDataHsh_ref", 10, "\n");#->1595
				next;
			}
			my $geneName = $DEResultDataHsh_ref->{$comparison}{$geneID}{'geneName'};
			@outputAry = ($geneID, $geneName);
			push @outputAry, @{$DEResultDataHsh_ref->{$comparison}{$geneID}{'fullDEInfo'}};
			foreach my $enrichmentInfoType (sort keys %{$geneBasedEnrichmentHsh_ref->{$comparison}}) {
				my $enrichmentStr = '--';
				$enrichmentStr = join ",", @{$geneBasedEnrichmentHsh_ref->{$comparison}{$enrichmentInfoType}{$geneID}} if exists $geneBasedEnrichmentHsh_ref->{$comparison}{$enrichmentInfoType}{$geneID};
				push @outputAry, $enrichmentStr;
			}
			print GENEBASEDSUM join "", (join "\t", @outputAry), "\n";
		}
		close GENEBASEDSUM;
	}
	
	return ($geneBasedSummaryPathHsh_ref);
}
sub printGeneListGAF {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: reportStatus|1595
#	appearInSub: runMap2SlimInAllComparison|1970
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_map2Slim|212
#	input: $geneGAFOriginalLineHsh_ref, $geneListGAFPath, $sgnfcntDEGeneNameHsh_ref
#	output: none
#	toCall: &printGeneListGAF($geneGAFOriginalLineHsh_ref, $sgnfcntDEGeneNameHsh_ref, $geneListGAFPath);
#	calledInLine: 914, 2008
#....................................................................................................................................................#

	#--&printGeneListGAF($geneGAFOriginalLineHsh_ref, $geneListHsh_ref, $geneListGAFPath);#->902

	my ($geneGAFOriginalLineHsh_ref, $sgnfcntDEGeneNameHsh_ref, $geneListGAFPath) = @_;
	
	&reportStatus("Printing geneListGAFPath", 10, "\n");#->1595

	open (GENELISTGAF, ">$geneListGAFPath");
	foreach my $geneID (sort {$a cmp $b} keys %{$sgnfcntDEGeneNameHsh_ref}) {
		if ($geneGAFOriginalLineHsh_ref->{$geneID}) {
			foreach my $theLine (@{$geneGAFOriginalLineHsh_ref->{$geneID}}) {
				print GENELISTGAF $theLine."\n" 
			}
		}
	}
	close GENELISTGAF;
}
sub readGODefObo {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: >none
#	appearInSub: runAndParseMap2Slim|1616
#	primaryAppearInSection: 4_processInputData|161
#	secondaryAppearInSection: >none
#	input: $GODefOboPath
#	output: $GODefOboHsh_ref, $GOSubsetHsh_ref
#	toCall: my ($GODefOboHsh_ref, $GOSubsetHsh_ref) = &readGODefObo($GODefOboPath);
#	calledInLine: 167, 1671
#....................................................................................................................................................#

	my ($GODefOboPath) = @_;
	
	open (GODEFOBO, "$GODefOboPath");
	
	my $GOID = "";
	my $object = "";
	my $GODefOboHsh_ref = {};
	my $GOSubsetHsh_ref = {};
	
	while (my $theLine = <GODEFOBO>) {
		$object = $1 if ($theLine =~ m/^\[(.+)\]$/);
		if ($object eq "Term") {
			$GOID = $1 if ($theLine =~ m/^id: (.+)$/);
			$GODefOboHsh_ref->{$GOID}{"name"} = $1 if ($theLine =~ m/^name: (.+)$/);
			$GODefOboHsh_ref->{$GOID}{"def"} = $1 if ($theLine =~ m/^def: (.+)$/);
			$GODefOboHsh_ref->{$GOID}{"namespace"} = $1 if ($theLine =~ m/^namespace: (.+)$/);
			$GOSubsetHsh_ref->{$1}{$GOID}++ if ($theLine =~ m/^subset: (.+)$/);
		}
	}
	close GODEFOBO;
	
	my $GONum = keys %{$GODefOboHsh_ref};
	my $subsetNum = keys %{$GOSubsetHsh_ref};
	
	return ($GODefOboHsh_ref, $GOSubsetHsh_ref);
}
sub readGeneDEResultInComparison {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|161
#	secondaryAppearInSection: >none
#	input: $DEComparisonInfoHsh_ref, $DEResultColumnIndexHsh_ref, $FCCutoffAry_ref, $geneInfoHsh_ref
#	output: $DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref, $sgnfcntDEGeneListHsh_ref
#	toCall: my ($sgnfcntDEGeneListHsh_ref, $DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref) = &readGeneDEResultInComparison($DEComparisonInfoHsh_ref, $FCCutoffAry_ref, $DEResultColumnIndexHsh_ref, $geneInfoHsh_ref);
#	calledInLine: 171
#....................................................................................................................................................#

	my ($DEComparisonInfoHsh_ref, $FCCutoffAry_ref , $DEResultColumnIndexHsh_ref, $geneInfoHsh_ref) = @_;

	my (%DEResultDataHsh, %sgnfcntDEGeneListHsh);
	
	my $fullDEInfoHeaderAry_ref = [];
	my $DEResultDataHsh_ref = {};
	my $sgnfcntDEGeneListHsh_ref = {};
	
	foreach my $comparison (sort {$a cmp $b} keys %{$DEComparisonInfoHsh_ref}) {
		&reportStatus("Reading DE results from $comparison", 10, "\n");#->1595

		my $DEResultPath = $DEComparisonInfoHsh_ref->{$comparison}{'DEResultPath'};
		my $refSample = $DEComparisonInfoHsh_ref->{$comparison}{'refSample'};
		my $qrySample = $DEComparisonInfoHsh_ref->{$comparison}{'qrySample'};
		
		open (DERESULT, "<", "$DEResultPath");
		my $header = <DERESULT>;
		while (my $theLine = <DERESULT>) {
			chomp $theLine;
			my @theLineSplt = split /\t/, $theLine;

			#---get gene ID and gene name
			my $geneID = $theLineSplt[0]; #---geneID has to be the first column
			my $geneName = 'unknown';
			$geneName = $geneInfoHsh_ref->{$geneID}{'geneName'} if $geneInfoHsh_ref->{$geneID}{'geneName'} ;
			$DEResultDataHsh_ref->{$comparison}{$geneID}{'geneName'} = $geneName;
			$DEResultDataHsh_ref->{$comparison}{$geneID}{'fullDEInfo'} = [];
			$fullDEInfoHeaderAry_ref = [];

			foreach my $item (sort keys %{$DEResultColumnIndexHsh_ref}) {
				my $value = $theLineSplt[$DEResultColumnIndexHsh_ref->{$item}];
				if ($value =~ m/^-?\d+\.?\d*$/ and $value ne '0') {#----http://docstore.mik.ua/orelly/perl4/cook/ch02_02.htm
					$value = sprintf "%.5f", $value;
				} else {
					if ($value eq '-inf') {
						$value = -9999;
					} elsif ($value eq 'inf') {
						$value = 9999;
					}
				}
				$DEResultDataHsh_ref->{$comparison}{$geneID}{$item} = $value;
				push @{$fullDEInfoHeaderAry_ref}, $item;
				push @{$DEResultDataHsh_ref->{$comparison}{$geneID}{'fullDEInfo'}}, $value;
			}

			foreach my $FCCutoff (@{$FCCutoffAry_ref}) {
				if (($DEResultDataHsh_ref->{$comparison}{$geneID}{'significantDiffExp'} ne '0') and (($DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_linearFC'} >= $FCCutoff) or ($DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_linearFC'} <= -$FCCutoff))) {
					$sgnfcntDEGeneListHsh_ref->{$comparison}{"FC$FCCutoff"}{$geneID} = $DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_linearFC'};
					print DEBUGLOG $FCCutoff."\t".$geneID."\t".$comparison."\t".$DEResultDataHsh_ref->{$comparison}{$geneID}{'significantDiffExp'}."\n";
				}
			}
		}
		close DERESULT;

		&reportStatus("$comparison results stored", 10, "\n");#->1595
	}
	
	return ($sgnfcntDEGeneListHsh_ref, $DEResultDataHsh_ref, $fullDEInfoHeaderAry_ref);
}
sub readGeneDEResultPathFile {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|161
#	secondaryAppearInSection: >none
#	input: $geneDEResultPathFile
#	output: $DEComparisonInfoHsh_ref
#	toCall: my ($DEComparisonInfoHsh_ref) = &readGeneDEResultPathFile($geneDEResultPathFile);
#	calledInLine: 170
#....................................................................................................................................................#
	
	my ($geneDEResultPathFile) = @_;
	
	my $DEComparisonInfoHsh_ref = {};
	
	open (DERESULT, "$geneDEResultPathFile");
	while (my $theLine = <DERESULT>) {
		chomp $theLine;
		my ($comparison, $refSample, $qrySample, $DEResultPath) = split /\t/, $theLine;
		$DEComparisonInfoHsh_ref->{$comparison}{'DEResultPath'} = $DEResultPath;
		die "Cant read $DEResultPath" if not -s $DEResultPath;
		&reportStatus("$comparison results checked", 10, "\n");#->1595
		$DEComparisonInfoHsh_ref->{$comparison}{'refSample'} = $refSample;
		$DEComparisonInfoHsh_ref->{$comparison}{'qrySample'} = $qrySample;
	}
	close DERESULT;
	
	return ($DEComparisonInfoHsh_ref);
}
sub readGeneGAF {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: >none
#	appearInSub: runAndParseMap2Slim|1616
#	primaryAppearInSection: 4_processInputData|161
#	secondaryAppearInSection: >none
#	input: $geneGAFPath
#	output: $geneGAFInfoHsh_ref, $geneGAFOriginalLineHsh_ref, $geneGOProcessHsh_ref
#	toCall: my ($geneGAFInfoHsh_ref, $geneGAFOriginalLineHsh_ref, $geneGOProcessHsh_ref) = &readGeneGAF($geneGAFPath);
#	calledInLine: 168, 1683
#....................................................................................................................................................#

	my ($geneGAFPath) = @_;
	
	my $geneGAFInfoHsh_ref = {};
	my $geneGAFOriginalLineHsh_ref = {};
	my $geneGOProcessHsh_ref = {};
	
	open (GENEGAF, "$geneGAFPath");
	while (my $theLine = <GENEGAF>) {
		chomp $theLine;
		my @theLineSplt = split /\t/, $theLine;
		my $geneID = $theLineSplt[1];
		my $GOID = $theLineSplt[4];
		my $evidence = $theLineSplt[6];
		my $GOProcess = $theLineSplt[8];
		
		$geneGOProcessHsh_ref->{$GOProcess}{$geneID}{$GOID}++;
		$geneGAFInfoHsh_ref->{$geneID}{$GOID} = $evidence;
		push @{$geneGAFOriginalLineHsh_ref->{$geneID}}, $theLine;
	}
	close GENEGAF;
	
	return ($geneGAFInfoHsh_ref, $geneGAFOriginalLineHsh_ref, $geneGOProcessHsh_ref);
}
sub readGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: getTextInfo
#	dependOnSub: reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|161
#	secondaryAppearInSection: >none
#	input: $geneInfoPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGeneInfo($geneInfoPath);
#	calledInLine: 169
#....................................................................................................................................................#

	my ($geneInfoPath)= @_;
	
	my $geneInfoHsh_ref = {};
	
	open (GENEINFO, "$geneInfoPath");
	while (my $theLine = <GENEINFO>) {
		chomp $theLine;
		my ($geneID, $geneName) = split /\t/, $theLine;
		$geneInfoHsh_ref->{$geneID}{'geneName'} = $geneName;
	}
	close GENEINFO;
	
	my $geneNum = keys %{$geneInfoHsh_ref};
	
	&reportStatus("$geneNum gene info read", 10, "\n");#->1595
	
	return ($geneInfoHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|66
#	secondaryAppearInSection: >none
#	input: none
#	output: $geneDEResultPathFile, $geneInfoPath, $outDir
#	toCall: my ($geneInfoPath, $geneDEResultPathFile, $outDir) = &readParameters();
#	calledInLine: 73
#....................................................................................................................................................#
	
	my ($geneInfoPath, $geneDEResultPathFile, $outDir);
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d_%02d_%02d_%02dhr%02dmn", $year+1900, $mon+1,$mday,$hour,$min;
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/$runTime/";


	GetOptions 	("geneInfoPath=s"  => \$geneInfoPath,
				 "geneDEResultPathFile=s"  => \$geneDEResultPathFile,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");


	#---check file
	foreach my $fileToCheck ($geneInfoPath, $geneDEResultPathFile) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}
	
	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -p -m 777 $outDir/";
	
	return ($geneInfoPath, $geneDEResultPathFile, $outDir);
}
sub refineAndParseFUNCHyperResults {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: generateGeneListHtmlLinkForSVG|536, generateRamiGOTreesvg|750, generate_random_string|790, reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 6_FUNCHyper|186
#	secondaryAppearInSection: >none
#	input: $DEComparisonInfoHsh_ref, $DEResultDataHsh_ref, $FDRCutOffOut, $FUNCHyperDir, $FUNCHyperRunInfoHsh_ref, $geneBasedEnrichmentHsh_ref, $geneInfoHsh_ref, $geneURLPrefix, $refinePVal, $rglrGODefOboHsh_ref, $sgnfcntDEGeneListHsh_ref
#	output: $FUNCHyperResultPathHsh_ref
#	toCall: my ($FUNCHyperResultPathHsh_ref) = &refineAndParseFUNCHyperResults($FUNCHyperRunInfoHsh_ref, $refinePVal, $FDRCutOffOut, $geneInfoHsh_ref, $DEResultDataHsh_ref, $sgnfcntDEGeneListHsh_ref, $FUNCHyperDir, $geneBasedEnrichmentHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEComparisonInfoHsh_ref);
#	calledInLine: 193
#....................................................................................................................................................#

	my ($FUNCHyperRunInfoHsh_ref, $refinePVal, $FDRCutOffOut, $geneInfoHsh_ref, $DEResultDataHsh_ref, $sgnfcntDEGeneListHsh_ref, $FUNCHyperDir, $geneBasedEnrichmentHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEComparisonInfoHsh_ref) = @_;
	
	#$FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{'runFUNCHyperCmd'} = $runFUNCHyperCmd;

	my $refinementLeftPVal = my $refinementRightPVal = $refinePVal;

	my $FUNCHyperResultDataHsh_ref = {};
	my $allGenesAtNodeHsh_ref = {};
	my $FUNCHyperResultPathHsh_ref = {};

	#-----define color for RamiGO
	my $colorHsh_ref = {};

	$colorHsh_ref->{'all'}{'+'} = "cyan4";
	$colorHsh_ref->{'all'}{'-'} = "cyan1";
	$colorHsh_ref->{'up'}{'+'} = "lawngreen";
	$colorHsh_ref->{'up'}{'-'} = "lightgreen";
	$colorHsh_ref->{'dn'}{'+'} = "red1";
	$colorHsh_ref->{'dn'}{'-'} = "pink";

	system "mkdir -m 777 $FUNCHyperDir/summary";

	foreach my $comparison (keys %{$FUNCHyperRunInfoHsh_ref}) {

		my $refSample = $DEComparisonInfoHsh_ref->{$comparison}{'refSample'};
		my $qrySample = $DEComparisonInfoHsh_ref->{$comparison}{'qrySample'};

		foreach my $cutoff (keys %{$FUNCHyperRunInfoHsh_ref->{$comparison}}) {
			
			my $GOColorRamiGOHsh_ref = {};
			
			my $mergeAllUpDnGOLinePath = "$FUNCHyperDir/summary/$comparison.$cutoff.merged.OneGOPerLine.tsv";
			my $mergeAllUpDnGeneLinePath = "$FUNCHyperDir/summary/$comparison.$cutoff.merged.OneGenePerLine.tsv";
			open (MERGEGO, '>', $mergeAllUpDnGOLinePath);
			open (MERGEGENE, '>', $mergeAllUpDnGeneLinePath);
			
			my $geneIDInGOGeneSetHsh_ref = {};
			
			foreach my $allUpDn (qw /up dn all/) {
				
				my $FUNCHyperOutDir = $FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{$allUpDn}{'FUNCHyperOutDir'};
				my $funcGOToGeneOutputPath = $FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{$allUpDn}{'funcGOToGeneOutputPath'};
				my $outputFUNCGOPath = "$FUNCHyperOutDir/groups.txt";
				&reportStatus("Parsing and refining FUNC results of $allUpDn $comparison at $cutoff", 10, "\n");#->1595
				
				#---get geneIN under GOID
				open (FUNCGOTOGENE, '<', $funcGOToGeneOutputPath);
				while (<FUNCGOTOGENE>) {
					chomp $_;
					my @splt = split /\t+/;
					my $geneID = $splt[0];
					foreach my $i (1..$#splt) {
						my $node_id = $splt[$i];
						next if $node_id eq 'all';
						#---get the sig_genes_in_node rather than total_genes_in_node;
						if (exists $sgnfcntDEGeneListHsh_ref->{$comparison}{$cutoff}{$geneID}) {
							my $DESeq1_linearFC = $sgnfcntDEGeneListHsh_ref->{$comparison}{$cutoff}{$geneID};
							$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'geneID'}{$geneID}++ 
							if (($DESeq1_linearFC > 0 and $allUpDn eq 'up')
							or	($DESeq1_linearFC < 0 and $allUpDn eq 'dn')
							or	($allUpDn eq 'all'));
						}
						$allGenesAtNodeHsh_ref->{$node_id}{'geneID'}{$geneID}++;
					}
				}
				close FUNCGOTOGENE;
				
				#---read the main results
				open (OUTPUTFUNCGO, '<', $outputFUNCGOPath);
				<OUTPUTFUNCGO>;
				while (<OUTPUTFUNCGO>) {
					chomp $_;
					my ($root_node_name, $node_name, $node_id, $total_genes_in_root, $total_genes_in_node, $sig_genes_in_root, $sig_genes_in_node, $raw_p_under, $raw_p_over, $FWER_under, $FWER_over, $FDR_under, $FDR_over) = split /\t+/;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'root_node_name'} = $root_node_name;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'node_name'} = $node_name;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'total_genes_in_root'} = $total_genes_in_root;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'total_genes_in_node'} = $total_genes_in_node;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'sig_genes_in_root'} = $sig_genes_in_root;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'sig_genes_in_node'} = $sig_genes_in_node;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'raw_p_under'} = $raw_p_under;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'raw_p_over'} = $raw_p_over;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'FWER_under'} = $FWER_under;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'FWER_over'} = $FWER_over;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'FDR_under'} = $FDR_under;
					$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'FDR_over'} = $FDR_over;
				}
				close OUTPUTFUNCGO;
				
				#---Run the refinement and parse them
				my ($shFile) = glob "$FUNCHyperOutDir/*.sh";
				my $refinementStatPath = "$FUNCHyperOutDir/refinement.stats.txt";
				system "$shFile $refinementLeftPVal $refinementRightPVal 1>$refinementStatPath 2>$refinementStatPath";
				my @refinedGOTerms = glob "$FUNCHyperOutDir/*$refinementLeftPVal\_$refinementRightPVal.txt";
				foreach my $refinedGOTermFiles (@refinedGOTerms) {
					open (REFINEGO, '<', $refinedGOTermFiles);
					<REFINEGO>;
					while (<REFINEGO>) {
						chomp $_;
						#cellular_component	external encapsulating structure	GO:0030312	-	0.853473	1	-1	-1	
						my ($root_node_name, $node_name, $node_id, $refine_sign, $raw_p_under, $raw_p_over, $refined_p_under, $refined_p_over) = split /\t+/;
						next unless exists $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id};
						$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'refine_sign'} = $refine_sign;
						$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'refined_p_under'} = $refined_p_under;
						$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'refined_p_over'} = $refined_p_over;
					}
					close REFINEGO;
				}
				
				{#---Parse all results
	
					my @headerAry = qw /root_node_name node_name total_genes_in_root total_genes_in_node sig_genes_in_root sig_genes_in_node FDR_over refined_p_over refine_sign/;

					my $random_string = &generate_random_string(5);#->790

					my $geneSetForGESAFile = "$FUNCHyperOutDir/GeneSetAtNodeForGESA.gmt";
					my $parsedOutputOneGOPerLineFile = "$FUNCHyperOutDir/pasred.overRep.OneGOPerLine.$random_string.tsv";
					open (PARSEDOUTGOLINE, '>', $parsedOutputOneGOPerLineFile);
					open (GENESETGSEA, '>', $geneSetForGESAFile);
					my @outputGOLineAry = ('node_id', 'allUpDn');
					push @outputGOLineAry, @headerAry;
					push @outputGOLineAry, 'gene';
					print PARSEDOUTGOLINE join "", ((join "\t",@outputGOLineAry), "\n");
					print MERGEGO join "", ((join "\t",@outputGOLineAry), "\n") if $allUpDn eq 'up'; #---make sure only print once
	
					my $parsedOutputOneGenePerLineFile = "$FUNCHyperOutDir/pasred.overRep.OneGenePerLine.$random_string.tsv";
					open (PARSEDOUTGENELINE, '>', $parsedOutputOneGenePerLineFile);
					my @outputGeneLineAry = ('node_id', 'allUpDn');
					push @outputGeneLineAry, @headerAry;
					push @outputGeneLineAry, (qw/geneID geneName DESeq1_linearFC signDE/);
					print PARSEDOUTGENELINE join "", ((join "\t",@outputGeneLineAry), "\n");
					print MERGEGENE join "", ((join "\t",@outputGeneLineAry), "\n") if $allUpDn eq 'up'; #---make sure only print once
					
					foreach my $node_id (sort {$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$a}{'FDR_over'} <=> $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$b}{'FDR_over'}} keys %{$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}}) {

						my $nodeName = $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'node_name'};

						#----not print the 'all' enrichment results if up or dn are enriched;
						my $printAllMergeResult = 1;
						if (($allUpDn eq 'all') and (exists $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{'dn'}{$node_id} or exists $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{'up'}{$node_id})) {
							$printAllMergeResult = 0;
						}
						
						#----ad hoc code to produce node GO geneset for GSEA
						my $geneNumInGeneSet =  keys %{$allGenesAtNodeHsh_ref->{$node_id}{'geneID'}};
						if ($geneNumInGeneSet >= 5) {
							my $geneSetTitle = "node_".$nodeName."_"."N=$geneNumInGeneSet"."_".$node_id;
							$geneSetTitle =~ s/\W+/\_/g;
							my @outputGeneSetGSEALineAry = ($geneSetTitle, $nodeName);
							foreach my $geneID (sort keys %{$allGenesAtNodeHsh_ref->{$node_id}{'geneID'}}) {
								push @outputGeneSetGSEALineAry, $geneID;
							}
							print GENESETGSEA join "", ((join "\t",@outputGeneSetGSEALineAry), "\n");
						}
	
						#----check FDR cutoff
						my $FDR_over = $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'FDR_over'};
						next if ($FDR_over < 0 or $FDR_over > $FDRCutOffOut);
						@outputGOLineAry = ($node_id, $allUpDn);
						foreach my $header (@headerAry) {
							push @outputGOLineAry, $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{$header};
						}
						
						foreach my $geneID (sort keys %{$allGenesAtNodeHsh_ref->{$node_id}{'geneID'}}) {
							my $geneName = 'unknown';
							$geneName = $geneInfoHsh_ref->{$geneID}{'geneName'} if $geneInfoHsh_ref->{$geneID}{'geneName'};
							@outputGeneLineAry = @outputGOLineAry;
							my $DESeq1_linearFC = $DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_linearFC'};
							my $signDE = 'No';
							if (exists $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'geneID'}{$geneID}) {
								$signDE = 'Yes';
								push @{$geneBasedEnrichmentHsh_ref->{$comparison}{"FuncHyper_all_".$cutoff}{$geneID}}, "$node_id:$nodeName";
								push @{$geneBasedEnrichmentHsh_ref->{$comparison}{"FuncHyper_refined".$cutoff}{$geneID}}, "$node_id:$nodeName" if $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'refine_sign'} eq '+';
							}
							
							$geneIDInGOGeneSetHsh_ref->{$node_id}{'gene'}{$geneID} = $signDE;

							push @outputGeneLineAry, ($geneID, $geneName, $DESeq1_linearFC, $signDE);
	
							print PARSEDOUTGENELINE join "", ((join "\t",@outputGeneLineAry), "\n");
							print MERGEGENE join "", ((join "\t",@outputGeneLineAry), "\n") if $printAllMergeResult;
						}
	
						my @tmpGeneIDAry = ();
						foreach my $geneID (sort keys %{$FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'geneID'}}) {
							push @tmpGeneIDAry, $geneID;
						}
						
						push @outputGOLineAry, join ";", sort @tmpGeneIDAry;
						print PARSEDOUTGOLINE join "", ((join "\t",@outputGOLineAry), "\n");
						
						if ($printAllMergeResult) {
							print MERGEGO join "", ((join "\t",@outputGOLineAry), "\n");
							my $GOID = $node_id;
							my $refine_sign = $FUNCHyperResultDataHsh_ref->{$comparison}{$cutoff}{$allUpDn}{$node_id}{'refine_sign'};
							my $color = $colorHsh_ref->{$allUpDn}{$refine_sign};
							$geneIDInGOGeneSetHsh_ref->{$GOID}{'upOrDn'} = $allUpDn;

							#----$GOColorRamiGOHsh_ref->{$GOID} already exists mean it appears in both up and down;
							if (exists $GOColorRamiGOHsh_ref->{$GOID}) {
								$color = $colorHsh_ref->{'all'}{$refine_sign};
								$geneIDInGOGeneSetHsh_ref->{$GOID}{'upOrDn'} = 'upAndDn';
							}

							$GOColorRamiGOHsh_ref->{$GOID} = $color;
						}
					}
					
					close PARSEDOUTGOLINE;
					close PARSEDOUTGENELINE;
					close GENESETGSEA;
				}
			}
			close MERGEGO;
			close MERGEGENE;
			
			my $RamiGOOutDir = "$FUNCHyperDir/summary/";
			my $svgName = "$comparison.$cutoff.merged";
			&generateRamiGOTreesvg($GOColorRamiGOHsh_ref, $RamiGOOutDir, $svgName);#->750
			my $svgPath = "$RamiGOOutDir"."$svgName.svg";
			my $svgDir = $RamiGOOutDir;
			my $analysisType = "FUNCHyper GO enrichment";
			&generateGeneListHtmlLinkForSVG($geneIDInGOGeneSetHsh_ref, $svgPath, $svgName, $svgDir, $comparison, $refSample, $qrySample, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEResultDataHsh_ref, $analysisType);#->536
			
			$FUNCHyperResultPathHsh_ref->{$cutoff}{$comparison}{'GOTsvPath'} = $mergeAllUpDnGOLinePath;
			$FUNCHyperResultPathHsh_ref->{$cutoff}{$comparison}{'geneTsvPath'} = $mergeAllUpDnGeneLinePath;
			$FUNCHyperResultPathHsh_ref->{$cutoff}{$comparison}{'svgPath'} = $svgPath;
			
			#---make sure the title "FuncHyper_".$cutoff exists in key
			$geneBasedEnrichmentHsh_ref->{$comparison}{"FuncHyper_all_".$cutoff} = {} if not exists $geneBasedEnrichmentHsh_ref->{$comparison}{"FuncHyper_all_".$cutoff};
			$geneBasedEnrichmentHsh_ref->{$comparison}{"FuncHyper_refined".$cutoff} = {} if not exists $geneBasedEnrichmentHsh_ref->{$comparison}{"FuncHyper_refined".$cutoff};
			
		}
	}
	
	return $FUNCHyperResultPathHsh_ref;
}
sub refineAndParseGSEAResults {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: generateGeneListHtmlLinkForSVG|536, generateRamiGOTreesvg|750, reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 7_GSEA|199
#	secondaryAppearInSection: >none
#	input: $DEComparisonInfoHsh_ref, $DEResultDataHsh_ref, $GSEARunInfoHsh_ref, $geneBasedEnrichmentHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref
#	output: $GSEAGOSVGPathHsh_ref, $GSEAResultDataHsh_ref, $GSEAResultPathHsh_ref
#	toCall: my ($GSEAResultDataHsh_ref, $GSEAResultPathHsh_ref, $GSEAGOSVGPathHsh_ref) = &refineAndParseGSEAResults($GSEARunInfoHsh_ref, $geneBasedEnrichmentHsh_ref, $DEResultDataHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEComparisonInfoHsh_ref);
#	calledInLine: 206
#....................................................................................................................................................#

	my ($GSEARunInfoHsh_ref, $geneBasedEnrichmentHsh_ref, $DEResultDataHsh_ref, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEComparisonInfoHsh_ref) = @_;
	
	my $colorHsh_ref = {};
	
	$colorHsh_ref->{'up'}{'hi'} = "lawngreen";
	$colorHsh_ref->{'up'}{'low'} = "lightgreen";
	$colorHsh_ref->{'dn'}{'hi'} = "red1";
	$colorHsh_ref->{'dn'}{'low'} = "pink";

	my $GSEAResultDataHsh_ref = {};
	my $GSEAResultPathHsh_ref = {};
	my $GSEAGOSVGPathHsh_ref = {};
	
	foreach my $comparison (sort keys %{$GSEARunInfoHsh_ref}) {

		my $refSample = $DEComparisonInfoHsh_ref->{$comparison}{'refSample'};
		my $qrySample = $DEComparisonInfoHsh_ref->{$comparison}{'qrySample'};
		
		foreach my $geneSetType ('GO', 'KEGG', 'FAM', 'INTERPRO') {

			&reportStatus("Parsing GSEA of $comparison", 10, "\n");#->1595

			my $GSEAOutDir = $GSEARunInfoHsh_ref->{$comparison}{$geneSetType}{'GSEAOutDir'};
			my ($GseaPrerankDir) = glob "$GSEAOutDir/*GseaPreranked*";
			my ($dnGeneSetXls) = glob "$GseaPrerankDir/gsea_report_for_na_neg*.xls";
			my ($upGeneSetXls) = glob "$GseaPrerankDir/gsea_report_for_na_pos*.xls";
			my ($indexHtml) = glob "$GseaPrerankDir/index.html";
			my @allHtmlAry = glob "$GseaPrerankDir/*.html";
			
			$GSEAResultPathHsh_ref->{$geneSetType}{$comparison}{'up'} = $upGeneSetXls;
			$GSEAResultPathHsh_ref->{$geneSetType}{$comparison}{'dn'} = $dnGeneSetXls;
			$GSEAResultPathHsh_ref->{$geneSetType}{$comparison}{'idx'} = $indexHtml;
		
			my $GSEAGOGeneSetHsh_ref = {};
		
			#-----change the lables of the HTMLs 
			foreach my $htmlFilePath (@allHtmlAry) {
				open (HTMLIN, "<", $htmlFilePath);
				my @bufferAry = <HTMLIN>;
				close HTMLIN;
				
				open (HTMLOUT, ">", $htmlFilePath);
				foreach (@bufferAry) {
					$_ =~ s/https\:\/\/www\.affymetrix\.com\/LinkServlet\?probeset\=/http\:\/\/amoebadb.org\/amoeba\/showRecord\.do\?name\=GeneRecordClasses\.GeneRecordClass&source_id\=/g;
					$_ =~ s/RANK METRIC SCORE/moderated log2 fold change/g;
					$_ =~ s/>na_pos</>$qrySample (comparing to $refSample)</g;
					$_ =~ s/>na_neg</>$refSample (comparing to $qrySample)</g;
					$_ =~ s/DESeq1_modLog2FC/$geneSetType for $comparison/g;
					
					if ($htmlFilePath =~ m/index.html/) {
						$_ =~ s/^<\/div><div><h4>Enrichment in phenotype: <b>na<\/b>/<\/div><div><h4>Upregulated in $qrySample/g;
						$_ =~ s/<\/ul><\/div><div><h4>Enrichment in phenotype: <b>na<\/b>/<\/div><div><h4>Upregulated in $refSample/g;
					}
					
					if ($htmlFilePath =~ m/neg_snapshot.html/) {
						$_ =~ s/Table: Snapshot of enrichment results/upregulated in $refSample of $geneSetType for $comparison/g;
						$_ =~ s/Table: Snapshot of enrichment results/upregulated in $refSample of $geneSetType for $comparison/g;
					}
					if ($htmlFilePath =~ m/pos_snapshot.html/) {
						$_ =~ s/Table: Snapshot of enrichment results/upregulated in $qrySample of $geneSetType for $comparison/g;
					}
					if ($htmlFilePath =~ m/gsea_report_for_na_neg_/) {
						$_ =~ s/Report for na_\w+ \d+ [GSEA]/upregulated in $refSample of $geneSetType for $comparison/g;
						$_ =~ s/Table: Gene sets enriched in phenotype <b>na<b>/upregulated in $refSample of $geneSetType for $comparison/g;
					}
					if ($htmlFilePath =~ m/gsea_report_for_na_pos_/) {
						$_ =~ s/Report for na_\w+ \d+ [GSEA]/upregulated in $qrySample of $geneSetType for $comparison/g;
						$_ =~ s/Table: Gene sets enriched in phenotype <b>na<b>/upregulated in $qrySample of $geneSetType for $comparison/g;
					}
					print HTMLOUT $_."\n";
				}
				close HTMLOUT;
			}

			foreach my $upOrDn ('up', 'dn') {
				my $resultPath = $GSEAResultPathHsh_ref->{$geneSetType}{$comparison}{$upOrDn};
				open (GSEARESULT, "<", $resultPath);
				<GSEARESULT>;
				while (<GSEARESULT>) {
					chomp $_;
					my ($geneSetName, $link, $details, $size, $ES, $NES, $p_val, $FDR, $FWER, $maxRnk, $leadingEdge) = split /\t/, $_;
					if ($geneSetName =~ m/GO.(\d+)/ and $geneSetType eq 'GO') {
						my $GOID = "GO:".$1;
						$GSEAGOGeneSetHsh_ref->{$upOrDn}{$GOID}{'NES'} = $NES;
						$GSEAGOGeneSetHsh_ref->{$upOrDn}{$GOID}{'geneSetName'} = $geneSetName;
						$GSEAGOGeneSetHsh_ref->{$upOrDn}{$GOID}{'size'} = $size;
						$GSEAGOGeneSetHsh_ref->{$upOrDn}{$GOID}{'FDR'} = $FDR;
					}
					$GSEAResultDataHsh_ref->{$comparison}{$geneSetType}{$geneSetName}{'upOrDn'} = $upOrDn;
					$GSEAResultDataHsh_ref->{$comparison}{$geneSetType}{$geneSetName}{'NES'} = $NES;
					$GSEAResultDataHsh_ref->{$comparison}{$geneSetType}{$geneSetName}{'size'} = $size;
					$GSEAResultDataHsh_ref->{$comparison}{$geneSetType}{$geneSetName}{'FDR'} = $FDR;
				}
				close GSEARESULT;
			}
		
			foreach my $FDRCutOff ((0.25, 0.1)) {
			
				my $geneIDInGOGeneSetHsh_ref = {};
				
				#---collect gene based data
				foreach my $geneSetName (keys %{$GSEAResultDataHsh_ref->{$comparison}{$geneSetType}}) {
					my $FDR = $GSEAResultDataHsh_ref->{$comparison}{$geneSetType}{$geneSetName}{'FDR'};
					next if  $FDR > $FDRCutOff;
					my ($theGeneSetXls) = glob "$GseaPrerankDir/$geneSetName.xls";
					open (GENESETXLS, "<", $theGeneSetXls);
					<GENESETXLS>;
					while (<GENESETXLS>) {
						chomp $_;
						my ($row, $geneID, $geneSymbol, $geneName, $rank, $score, $runningES, $coreEnrichment) = split /\t/, $_;
						if ($geneSetName =~ m/GO.(\d+)/ and $geneSetType eq 'GO') {
							my $GOID = "GO:".$1;
							$geneIDInGOGeneSetHsh_ref->{$GOID}{'gene'}{$geneID} = $coreEnrichment;
						}
						push @{$geneBasedEnrichmentHsh_ref->{$comparison}{"$geneSetType\_GSEA_all_FDR".$FDRCutOff}{$geneID}}, "$geneSetName";
						push @{$geneBasedEnrichmentHsh_ref->{$comparison}{"$geneSetType\_GSEA_core_FDR".$FDRCutOff}{$geneID}}, "$geneSetName" if $coreEnrichment eq 'Yes';
					}
					close GENESETXLS;
				}

				$geneBasedEnrichmentHsh_ref->{$comparison}{"$geneSetType\_GSEA_all_FDR".$FDRCutOff} = {} if not exists $geneBasedEnrichmentHsh_ref->{$comparison}{"$geneSetType\_GSEA_all_FDR".$FDRCutOff};
				$geneBasedEnrichmentHsh_ref->{$comparison}{"$geneSetType\_GSEA_core_FDR".$FDRCutOff} = {} if not exists $geneBasedEnrichmentHsh_ref->{$comparison}{"$geneSetType\_GSEA_core_FDR".$FDRCutOff};

				#----make GO graph
				if ($geneSetType eq 'GO') {
					my $GOColorRamiGOHsh_ref = {};
					foreach my $upOrDn (keys %{$GSEAGOGeneSetHsh_ref}) {
						foreach my $GOID (keys %{$GSEAGOGeneSetHsh_ref->{$upOrDn}}) {
							my $FDR = $GSEAGOGeneSetHsh_ref->{$upOrDn}{$GOID}{'FDR'};
							next if  $FDR > $FDRCutOff;
							my $FDRLevel = 'low';
							$FDRLevel = 'hi' if $FDR <= 0.1;
							my $color = $colorHsh_ref->{$upOrDn}{$FDRLevel};
							$GOColorRamiGOHsh_ref->{$GOID} = $color;
							$geneIDInGOGeneSetHsh_ref->{$GOID}{'upOrDn'} = $upOrDn;
						}
					}

					system ("mkdir -m 777 $GSEAOutDir/summary/ >/dev/null 2>&1");
					my $RamiGOOutDir = "$GSEAOutDir/summary/";
					my $svgName = "$comparison.FDR.$FDRCutOff";
					&generateRamiGOTreesvg($GOColorRamiGOHsh_ref, $RamiGOOutDir, $svgName);#->750
					my $svgPath = "$RamiGOOutDir"."$svgName.svg";
					my $svgDir = $RamiGOOutDir;
					my $analysisType = "GSEA GO Geneset";
					&generateGeneListHtmlLinkForSVG($geneIDInGOGeneSetHsh_ref, $svgPath, $svgName, $svgDir, $comparison, $refSample, $qrySample, $geneURLPrefix, $rglrGODefOboHsh_ref, $DEResultDataHsh_ref, $analysisType);#->536
					$GSEAGOSVGPathHsh_ref->{$comparison}{$FDRCutOff} = $svgPath;
				}
			}
		}
	}
	
	return $GSEAResultDataHsh_ref, $GSEAResultPathHsh_ref, $GSEAGOSVGPathHsh_ref;
	
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|280
#	appearInSub: generateRamiGOTreesvg|750, printGeneBasedSummary|851, printGeneListGAF|902, readGeneDEResultInComparison|970, readGeneDEResultPathFile|1042, readGeneInfo|1109, refineAndParseFUNCHyperResults|1178, refineAndParseGSEAResults|1427, runAndParseMap2Slim|1616, runFUNCHyperInAllComparison|1852, runGSEAInAllComparison|1917, runMap2SlimInAllComparison|1970
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|161, 6_FUNCHyper|186, 7_GSEA|199, 8_map2Slim|212, 9_outputData|225
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 767, 883, 918, 991, 1036, 1064, 1135, 1233, 1459, 1642, 1662, 1765, 1894, 1951, 2001
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->280

	return ();
}
sub runAndParseMap2Slim {
#....................................................................................................................................................#
#	subroutineCategory: GeneOntology
#	dependOnSub: readGODefObo|931, readGeneGAF|1073, reportStatus|1595
#	appearInSub: runMap2SlimInAllComparison|1970
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 8_map2Slim|212
#	input: $DEResultDataHsh_ref, $GOSlimOboPathAry_ref, $GOSlimOutDir, $comparison, $fullDEInfoHeaderAry_ref, $geneListGAFPath, $regularGOOboPath, $sgnfcntDEGeneNameHsh_ref
#	output: $GOBasedOneLinerGOSlimResultPath, $GOColorRamiGOHsh_ref, $geneBasedOneLinerGOSlimResultPath, $geneGOSlimOneToOneResultPath
#	toCall: my ($GOColorRamiGOHsh_ref, $geneBasedOneLinerGOSlimResultPath, $GOBasedOneLinerGOSlimResultPath, $geneGOSlimOneToOneResultPath) = &runAndParseMap2Slim($GOSlimOboPathAry_ref, $geneListGAFPath, $regularGOOboPath, $sgnfcntDEGeneNameHsh_ref, $GOSlimOutDir, $DEResultDataHsh_ref, $comparison, $fullDEInfoHeaderAry_ref);
#	calledInLine: 2009
#....................................................................................................................................................#

	my ($GOSlimOboPathAry_ref, $geneListGAFPath, $regularGOOboPath, $sgnfcntDEGeneNameHsh_ref, $GOSlimOutDir, $DEResultDataHsh_ref, $comparison, $fullDEInfoHeaderAry_ref) = @_;

	my %mergedGOSlimDefOboHsh;
	my %mergedGOSlimGAFInfoHsh;
	my %mergedGOSlimGOIDSourceHsh;
	
	my $geneNum = keys %{$sgnfcntDEGeneNameHsh_ref};

	#---issue GOSlim CMD in parallel 
	my $GOSlimOutPathHsh_ref = {};
	foreach my $GOSlimOboPath (@{$GOSlimOboPathAry_ref}) {

		my ($fileTag, undef, undef) = fileparse($GOSlimOboPath, qr/\.[^.]*/);

		&reportStatus("Running GOSLIM $fileTag in parallel", 10, "\n");#->1595

		system ("mkdir -p -m 777 $GOSlimOutDir/list/individual/$fileTag/");
		my $GOSlimOutPath = "$GOSlimOutDir/list/individual/$fileTag/map2Slim.gaf";
		$GOSlimOutPathHsh_ref->{$GOSlimOboPath} = $GOSlimOutPath;
		my $map2SlimCMD = "map2slim $GOSlimOboPath $regularGOOboPath $geneListGAFPath -o $GOSlimOutPath";
		my $errorLogPath = "$GOSlimOutDir/list/individual/$fileTag/GOSlim.error.log.txt";
		system (qq|$map2SlimCMD 1>$errorLogPath 2>$errorLogPath &|);
	}

	#---will not proceed if there's map2slim in ps
	my $grepCMD = "ps -ef | grep map2slim | grep -v grep";
	my $grepStr = "map2slim";
	my $sdout = $grepStr;
	while ($sdout =~ m/$grepStr/) {
		$sdout = `$grepCMD`;
		sleep (0.1);
	}

	#---will not proceed if there's map2slim in ps
	&reportStatus("Printing GOSlim results", 10, "\n");#->1595
	foreach my $GOSlimOboPath (@{$GOSlimOboPathAry_ref}) {

		my @filePathAry = split /\//, $GOSlimOboPath;
		my $fileTag = $filePathAry[-1];
		$fileTag =~ s/\.obo$//;

		#---read the GOSlim results
		my $GOSlimOutPath = $GOSlimOutPathHsh_ref->{$GOSlimOboPath};
		my ($GOSlimDefOboHsh_ref, $GOSlimSubsetHsh_ref)= &readGODefObo($GOSlimOboPath);#->931
		my %GOSlimDefOboHsh = %{$GOSlimDefOboHsh_ref};
		
		#----put the individual GOSlim Def into the merge GOSlim Def
		foreach my $GOID (keys %GOSlimDefOboHsh) {
			#---skipp the general GOIDs
			${$mergedGOSlimDefOboHsh{$GOID}}{"name"} = ${$GOSlimDefOboHsh{$GOID}}{"name"};
			${$mergedGOSlimDefOboHsh{$GOID}}{"def"} = ${$GOSlimDefOboHsh{$GOID}}{"def"};
			${$mergedGOSlimDefOboHsh{$GOID}}{"namespace"} = ${$GOSlimDefOboHsh{$GOID}}{"namespace"};
		}

		#----parse the GOSlim results
		my ($GOSlimGAFInfoHsh_ref, $GOSlimGAFOriginalLineHsh_ref, undef) = &readGeneGAF($GOSlimOutPath);#->1073
		my %GOSlimGAFInfoHsh = %{$GOSlimGAFInfoHsh_ref};
		
		my %GOBasedCountHsh;
		my $geneBasedOneLinerGOSlimResultPath = $GOSlimOutPath.".geneBased.oneLiner.n=$geneNum.log.tsv";
		my $geneGOSlimOneToOneResultPath = $GOSlimOutPath.".geneGOSlim.OneToOne.n=$geneNum.log.tsv";

		open (GENEBASEDONE, ">$geneBasedOneLinerGOSlimResultPath");
		print GENEBASEDONE join '', ((join "\t", ("geneID", "geneName", "GONum", "GOStr", @{$fullDEInfoHeaderAry_ref})), "\n");
		open (GENEGOONETOONE, ">$geneGOSlimOneToOneResultPath");
		
		print GENEGOONETOONE join '', ((join "\t", ("namespace", "geneID", "geneName", "GOID", "name", @{$fullDEInfoHeaderAry_ref}, "evidence", "def")), "\n");
		
		foreach my $geneID (sort {$a cmp $b} keys %{$sgnfcntDEGeneNameHsh_ref}) {
			my $geneName = $sgnfcntDEGeneNameHsh_ref->{$geneID};
			my @fullDEInfoAry = @{$DEResultDataHsh_ref->{$comparison}{$geneID}{'fullDEInfo'}};
			my @GOAry = ();
			if (exists $GOSlimGAFInfoHsh{$geneID}) {
				foreach my $GOID (sort {$a cmp $b} keys %{$GOSlimGAFInfoHsh{$geneID}}) {
					
					next if ($GOID eq 'GO:0005575' 
					or 	$GOID eq 'GO:0005623' 
					or 	$GOID eq 'GO:0008150' 
					or 	$GOID eq 'GO:0065007' 
					or 	$GOID eq 'GO:0009987' 
					or 	$GOID eq 'GO:0023052' 
					or 	$GOID eq 'GO:0003824'
					or 	$GOID eq 'GO:0003674');

					#--put the individual GOSlim result into the merge GOSlim result
					my $evidence = ${$GOSlimGAFInfoHsh{$geneID}}{$GOID};
					${$mergedGOSlimGAFInfoHsh{$geneID}}{$GOID} = $evidence;
					${$mergedGOSlimGOIDSourceHsh{$GOID}}{$fileTag}++;
					
					${$GOBasedCountHsh{$GOID}}{$geneID} = $geneName;
					my $name = ${$GOSlimDefOboHsh{$GOID}}{"name"};
					my $def = ${$GOSlimDefOboHsh{$GOID}}{"def"};
					my $namespace = ${$GOSlimDefOboHsh{$GOID}}{"namespace"};
					print GENEGOONETOONE join '', ((join "\t", ($namespace, $geneID, $geneName, $GOID, $name, @fullDEInfoAry, $evidence, $def)), "\n");
					push @GOAry, $GOID;
				}
			} else {
				my $namespace = my $GOID = my $name = my $def = my $evidence = '#null';
				print GENEGOONETOONE join '', ((join "\t", ($namespace, $geneID, $geneName, $GOID, $name, @fullDEInfoAry, $evidence, $def)), "\n");
			}
			push @GOAry, 'null' if @GOAry == 0;
			my $GONum = @GOAry;
			my $GOStr = join ";", @GOAry;
			print GENEBASEDONE join '', ((join "\t", ($geneID, $geneName, $GONum, $GOStr, @fullDEInfoAry)), "\n");
		}
		close GENEBASEDONE;
		close GENEGOONETOONE;
		
		my $GOBasedOneLinerGOSlimResultPath = $GOSlimOutPath.".GOBased.OneLiner.n=$geneNum.log.tsv";
		open (GOBASEDONE, ">$GOBasedOneLinerGOSlimResultPath");
		print GOBASEDONE join '', ((join "\t", ("namespace", "GOID", "name", "geneNum", "def", "geneStr")), "\n");
		foreach my $GOID (sort {$a cmp $b} keys %GOBasedCountHsh) {
			
			next if ($GOID eq 'GO:0005575' 
			or 	$GOID eq 'GO:0005623' 
			or 	$GOID eq 'GO:0008150' 
			or 	$GOID eq 'GO:0065007' 
			or 	$GOID eq 'GO:0009987' 
			or 	$GOID eq 'GO:0023052' 
			or 	$GOID eq 'GO:0003824'
			or 	$GOID eq 'GO:0003674');

			my $name = ${$GOSlimDefOboHsh{$GOID}}{"name"};
			my $def = ${$GOSlimDefOboHsh{$GOID}}{"def"};
			my $namespace = ${$GOSlimDefOboHsh{$GOID}}{"namespace"};
			my @geneIDAry;
			foreach my $geneID (sort {$a cmp $b} keys %{$GOBasedCountHsh{$GOID}}) {
				push @geneIDAry, $geneID;
			}
			my $geneNum = @geneIDAry;
			my $geneStr = join ";", @geneIDAry;
			print GOBASEDONE join '', ((join "\t", ($namespace, $GOID, $name, $geneNum, $def, $geneStr)), "\n");
		}
		close GOBASEDONE;
	}
	
	#----print the merged GOSlim results
	&reportStatus("Printing merged GOSlim results", 10, "\n");#->1595
	system ("mkdir -p -m 777 $GOSlimOutDir/list/merged/");
	my $GOSlimOutPath = "$GOSlimOutDir/list/merged/$comparison.merged";
	my %GOBasedCountHsh;
	my $geneBasedOneLinerGOSlimResultPath = $GOSlimOutPath.".geneBased.oneLiner.n=$geneNum.log.tsv";
	my $geneGOSlimOneToOneResultPath = $GOSlimOutPath.".geneGOSlim.OneToOne.n=$geneNum.log.tsv";
	open (GENEBASEDONE, ">$geneBasedOneLinerGOSlimResultPath");
	print GENEBASEDONE join '', ((join "\t", ("geneID", "geneName", "GONum", "GOStr", @{$fullDEInfoHeaderAry_ref})), "\n");
	open (GENEGOONETOONE, ">$geneGOSlimOneToOneResultPath");
	print GENEGOONETOONE join '', ((join "\t", ("namespace", "geneID", "geneName", "GOID", "name", "GOSlimSource", @{$fullDEInfoHeaderAry_ref}, "evidence", "def")), "\n");

	foreach my $geneID (sort {$a cmp $b} keys %{$sgnfcntDEGeneNameHsh_ref}) {
		my $geneName = $sgnfcntDEGeneNameHsh_ref->{$geneID};
		my @fullDEInfoAry = @{$DEResultDataHsh_ref->{$comparison}{$geneID}{'fullDEInfo'}};
		my @GOAry = ();
		if (exists $mergedGOSlimGAFInfoHsh{$geneID}) {
			foreach my $GOID (sort {$a cmp $b} keys %{$mergedGOSlimGAFInfoHsh{$geneID}}) {
				
				#--put the individual GOSlim result into the merge GOSlim result
				my @GOSlimSourceAry;
				foreach my $GOSlimSource (sort {${$mergedGOSlimGOIDSourceHsh{$GOID}}{$b} <=> ${$mergedGOSlimGOIDSourceHsh{$GOID}}{$a}} keys %{$mergedGOSlimGOIDSourceHsh{$GOID}}) {
					push @GOSlimSourceAry, $GOSlimSource."=".${$mergedGOSlimGOIDSourceHsh{$GOID}}{$GOSlimSource};
				}
				my $GOSlimSource = join ";", @GOSlimSourceAry;

				my $evidence = ${$mergedGOSlimGAFInfoHsh{$geneID}}{$GOID};
				${$GOBasedCountHsh{$GOID}}{$geneID} = $geneName;
				my $name = ${$mergedGOSlimDefOboHsh{$GOID}}{"name"};
				my $def = ${$mergedGOSlimDefOboHsh{$GOID}}{"def"};
				my $namespace = ${$mergedGOSlimDefOboHsh{$GOID}}{"namespace"};
				print GENEGOONETOONE join '', ((join "\t", ($namespace, $geneID, $geneName, $GOID, $name, $GOSlimSource, @fullDEInfoAry, $evidence, $def)), "\n");
				push @GOAry, $GOID;
			}
		} else {
			my $namespace = my $GOID = my $name = my $def = my $evidence = my $GOSlimSource = '#null';
			print GENEGOONETOONE join '', ((join "\t", ($namespace, $geneID, $geneName, $GOID, $name, $GOSlimSource, @fullDEInfoAry, $evidence, $def)), "\n");
		}
		
		my $GONum = @GOAry;
		my $GOStr = join ";", @GOAry;
		print GENEBASEDONE join '', ((join "\t", ($geneID, $geneName, $GONum, $GOStr, @fullDEInfoAry)), "\n");
	}
	close GENEBASEDONE;
	close GENEGOONETOONE;
	
	my $GOColorRamiGOHsh_ref = {};
	
	my $GOBasedOneLinerGOSlimResultPath = $GOSlimOutPath.".GOBased.OneLiner.n=$geneNum.log.tsv";
	open (GOBASEDONE, ">$GOBasedOneLinerGOSlimResultPath");
	print GOBASEDONE join '', ((join "\t", ("namespace", "GOSlimSource", "GOID", "name", "geneNum", "def", "geneStr")), "\n");
	foreach my $GOID (sort {$a cmp $b} keys %GOBasedCountHsh) {
		my @GOSlimSourceAry;
		my $FCDirHsh_ref = {};
		
		foreach my $GOSlimSource (sort {${$mergedGOSlimGOIDSourceHsh{$GOID}}{$b} <=> ${$mergedGOSlimGOIDSourceHsh{$GOID}}{$a}} keys %{$mergedGOSlimGOIDSourceHsh{$GOID}}) {
			push @GOSlimSourceAry, $GOSlimSource."=".${$mergedGOSlimGOIDSourceHsh{$GOID}}{$GOSlimSource};
		}
		my $GOSlimSource = join ";", @GOSlimSourceAry;

		my $name = ${$mergedGOSlimDefOboHsh{$GOID}}{"name"};
		my $def = ${$mergedGOSlimDefOboHsh{$GOID}}{"def"};
		my $namespace = ${$mergedGOSlimDefOboHsh{$GOID}}{"namespace"};
		my @geneIDAry;
		
		foreach my $geneID (sort {$a cmp $b} keys %{$GOBasedCountHsh{$GOID}}) {
			my $DESeq1_modLog2FC = $DEResultDataHsh_ref->{$comparison}{$geneID}{'DESeq1_modLog2FC'};
			$FCDirHsh_ref->{'up'}++ if $DESeq1_modLog2FC > 0;
			$FCDirHsh_ref->{'dn'}++ if $DESeq1_modLog2FC < 0;
			push @geneIDAry, $geneID;
		}
		
		my $RAMIGOColor;
		$RAMIGOColor = 'lightgreen' if $FCDirHsh_ref->{'up'} and not $FCDirHsh_ref->{'dn'};
		$RAMIGOColor = 'pink' if $FCDirHsh_ref->{'dn'} and not $FCDirHsh_ref->{'up'};
		$RAMIGOColor = 'gold' if $FCDirHsh_ref->{'dn'} and $FCDirHsh_ref->{'up'};
		$GOColorRamiGOHsh_ref->{$GOID} = $RAMIGOColor;
		
		my $geneNum = @geneIDAry;
		my $geneStr = join ";", @geneIDAry;
		print GOBASEDONE join '', ((join "\t", ($namespace, $GOSlimSource, $GOID, $name, $geneNum, $def, $geneStr)), "\n");
	}
	close GOBASEDONE;
	
	return $GOColorRamiGOHsh_ref, $geneBasedOneLinerGOSlimResultPath, $GOBasedOneLinerGOSlimResultPath, $geneGOSlimOneToOneResultPath;
	
}
sub runFUNCHyperInAllComparison {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: generateFUNCHyperInput|298, reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 6_FUNCHyper|186
#	secondaryAppearInSection: >none
#	input: $DEResultDataHsh_ref, $FPKMCutoff, $FUNCGOTermDBDir, $FUNCGeneToGOScript, $FUNCGraphPath, $FUNCHyperDir, $FUNCTermPath, $geneGAFInfoHsh_ref, $maxThread, $sgnfcntDEGeneListHsh_ref
#	output: $FUNCHyperRunInfoHsh_ref
#	toCall: my ($FUNCHyperRunInfoHsh_ref) = &runFUNCHyperInAllComparison($sgnfcntDEGeneListHsh_ref, $FUNCHyperDir, $geneGAFInfoHsh_ref, $FUNCGOTermDBDir, $FUNCTermPath, $FUNCGraphPath, $FUNCGeneToGOScript, $FPKMCutoff, $DEResultDataHsh_ref, $maxThread);
#	calledInLine: 192
#....................................................................................................................................................#

	my ($sgnfcntDEGeneListHsh_ref, $FUNCHyperDir, $geneGAFInfoHsh_ref, $FUNCGOTermDBDir, $FUNCTermPath, $FUNCGraphPath, $FUNCGeneToGOScript, $FPKMCutoff, $DEResultDataHsh_ref, $maxThread) = @_;

	my $threadInitiated = 0;
	my $totalComparisonNum = my $procComparisonNum = 0;
	foreach my $comparison (sort {$a cmp $b} keys %{$sgnfcntDEGeneListHsh_ref}) {
		foreach my $cutoff (sort {$a cmp $b} keys %{$sgnfcntDEGeneListHsh_ref->{$comparison}}) {
			foreach my $allUpDn (qw /all up dn/) {
				$totalComparisonNum++;
			}
		}
	}

	my $FUNCHyperRunInfoHsh_ref = {};

	foreach my $comparison (sort {$a cmp $b} keys %{$sgnfcntDEGeneListHsh_ref}) {
		foreach my $cutoff (sort {$a cmp $b} keys %{$sgnfcntDEGeneListHsh_ref->{$comparison}}) {
			foreach my $allUpDn (qw /all up dn/) {
				$procComparisonNum++;

				system ("mkdir -p -m 777 $FUNCHyperDir/$comparison/$cutoff/");
				my ($runFUNCHyperCmd, $FUNCHyperOutDir, $funcGOToGeneOutputPath, $funcGOToGeneCMD) =  &generateFUNCHyperInput($FUNCHyperDir, $sgnfcntDEGeneListHsh_ref, $geneGAFInfoHsh_ref, $comparison, $cutoff, $FUNCGOTermDBDir, $FUNCTermPath, $FUNCGraphPath, $FUNCGeneToGOScript, $allUpDn, $FPKMCutoff, $DEResultDataHsh_ref);#->298
				my $errorLogPath = "$FUNCHyperDir/$comparison/$cutoff/error.log.txt";
				
				$FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{$allUpDn}{'funcGOToGeneCMD'} = $funcGOToGeneCMD;
				$FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{$allUpDn}{'runFUNCHyperCmd'} = $runFUNCHyperCmd;
				$FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{$allUpDn}{'FUNCHyperOutDir'} = $FUNCHyperOutDir;
				$FUNCHyperRunInfoHsh_ref->{$comparison}{$cutoff}{$allUpDn}{'funcGOToGeneOutputPath'} = $funcGOToGeneOutputPath;
				
				system (qq|$runFUNCHyperCmd 1>$errorLogPath 2>$errorLogPath & $funcGOToGeneCMD &|);
	
				&reportStatus("Running runFUNCHyper for $allUpDn $comparison at $cutoff", 10, "\n");#->1595
	
				$threadInitiated++;

				if ($threadInitiated == $maxThread or $procComparisonNum == $totalComparisonNum) {#----hole the process while maxThread Running or last one
					$threadInitiated = 0;
					#---will not proceed if there's map2slim in ps
					my $grepCMD = "ps -ef | grep func_hyper | grep -v grep";
					my $grepStr = "func_hyper";
					my $sdout = $grepStr;
					while ($sdout =~ m/$grepStr/) {
						$sdout = `$grepCMD`;
						sleep (0.1);
					}
				}
			}
		}
	}

	return $FUNCHyperRunInfoHsh_ref;

}
sub runGSEAInAllComparison {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: generateGSEAInput|409, reportStatus|1595
#	appearInSub: >none
#	primaryAppearInSection: 7_GSEA|199
#	secondaryAppearInSection: >none
#	input: $DEResultDataHsh_ref, $FPKMCutoff, $GSEABinPath, $GSEADir, $geneSetPathHsh_ref, $maxGeneSetSize, $minGeneSetSize
#	output: $GSEARunInfoHsh_ref
#	toCall: my ($GSEARunInfoHsh_ref) = &runGSEAInAllComparison($DEResultDataHsh_ref, $GSEADir, $GSEABinPath, $maxGeneSetSize, $minGeneSetSize, $FPKMCutoff, $geneSetPathHsh_ref);
#	calledInLine: 205
#....................................................................................................................................................#

	my ($DEResultDataHsh_ref, $GSEADir, $GSEABinPath, $maxGeneSetSize, $minGeneSetSize, $FPKMCutoff, $geneSetPathHsh_ref, $maxThread) = @_;
	
	my $threadInitiated = 0;
	my $totalComparisonNum = keys %{$DEResultDataHsh_ref};
	my $procComparisonNum = 0;
	
	my $GSEARunInfoHsh_ref = {};
	
	foreach my $comparison (sort {$a cmp $b} keys %{$DEResultDataHsh_ref}) {
		$procComparisonNum++;
		system ("mkdir -p -m 777 $GSEADir/$comparison/");
		my ($GSEACmdAndDirHsh_ref) = &generateGSEAInput($GSEADir, $DEResultDataHsh_ref, $geneSetPathHsh_ref, $comparison, $GSEABinPath, $maxGeneSetSize, $minGeneSetSize, $FPKMCutoff);#->409
		my $errorLogPath = "$GSEADir/$comparison/error.log.txt";
		
		foreach my $geneSetType (keys %{$GSEACmdAndDirHsh_ref}) {
			$GSEARunInfoHsh_ref->{$comparison}{$geneSetType}{'runGSEACmd'} = $GSEACmdAndDirHsh_ref->{$geneSetType}{'runGSEACmd'};
			$GSEARunInfoHsh_ref->{$comparison}{$geneSetType}{'GSEAOutDir'} = $GSEACmdAndDirHsh_ref->{$geneSetType}{'GSEAOutDir'};
			system (qq|$GSEARunInfoHsh_ref->{$comparison}{$geneSetType}{'runGSEACmd'} 1>$errorLogPath 2>$errorLogPath &|);
			$threadInitiated++;
		}

		&reportStatus("Running GSEA for $comparison", 10, "\n");#->1595
		
		if (($threadInitiated >= $maxThread) or ($procComparisonNum == $totalComparisonNum)) {#----hole the process while maxThread Running or last one
			$threadInitiated = 0;
			#---will not proceed if there's map2slim in ps
			my $grepCMD = "ps -ef | grep xtools.gsea.GseaPreranked | grep -v grep";
			my $grepStr = "xtools.gsea.GseaPreranked";
			my $sdout = $grepStr;
			while ($sdout =~ m/$grepStr/) {
				$sdout = `$grepCMD`;
				sleep (0.1);
			}
		}
	}
	
	return $GSEARunInfoHsh_ref;
	
}
sub runMap2SlimInAllComparison {
#....................................................................................................................................................#
#	subroutineCategory: specific
#	dependOnSub: generateRamiGOTreesvg|750, printGeneListGAF|902, reportStatus|1595, runAndParseMap2Slim|1616
#	appearInSub: >none
#	primaryAppearInSection: 8_map2Slim|212
#	secondaryAppearInSection: >none
#	input: $DEResultDataHsh_ref, $GOSlimOboPathAry_ref, $fullDEInfoHeaderAry_ref, $geneGAFOriginalLineHsh_ref, $outDir, $regularGOOboPath, $sgnfcntDEGeneListHsh_ref
#	output: $GOSlimResultPathHsh_ref, $GOSlimRootDir
#	toCall: my ($GOSlimResultPathHsh_ref, $GOSlimRootDir) = &runMap2SlimInAllComparison($sgnfcntDEGeneListHsh_ref, $DEResultDataHsh_ref, $geneGAFOriginalLineHsh_ref, $GOSlimOboPathAry_ref, $regularGOOboPath, $fullDEInfoHeaderAry_ref, $outDir);
#	calledInLine: 218
#....................................................................................................................................................#

	my ($sgnfcntDEGeneListHsh_ref, $DEResultDataHsh_ref, $geneGAFOriginalLineHsh_ref, $GOSlimOboPathAry_ref, $regularGOOboPath, $fullDEInfoHeaderAry_ref, $outDir) = @_;

	my $GOSlimResultPathHsh_ref = {};
	my $GOSlimRootDir = "$outDir/GOSlim/";
	system ("mkdir -p -m 777 $outDir/overall/");
	open (OVERALLGENENUM, ">$outDir/overall/comparison.cutoff.geneNum.log.txt");
	foreach my $comparison (sort {$a cmp $b} keys %{$sgnfcntDEGeneListHsh_ref}) {
		system ("mkdir -p -m 777 $outDir/GOSlim/$comparison/");
		foreach my $cutoff (sort {$a cmp $b} keys %{$sgnfcntDEGeneListHsh_ref->{$comparison}}) {
			
			my $sgnfcntDEGeneNameHsh_ref = {}; 
			
			foreach my $geneID (keys %{$sgnfcntDEGeneListHsh_ref->{$comparison}{$cutoff}}) {
				$sgnfcntDEGeneNameHsh_ref->{$geneID} = $DEResultDataHsh_ref->{$comparison}{$geneID}{'geneName'};
			}
			
			my $geneNum = keys %{$sgnfcntDEGeneNameHsh_ref};
			my $GOSlimOutDir = "$GOSlimRootDir/$comparison/$cutoff.n=$geneNum/";
			system ("mkdir -p -m 777 $GOSlimOutDir");
			&reportStatus("GOSlim for $comparison at $cutoff n=$geneNum", 10, "\n");#->1595
			print OVERALLGENENUM "$comparison\t$cutoff\t$geneNum\n";
			if ($geneNum <= 15) {
				print "less than 15, skipped\n.";
				next;
			}
			my $geneListGAFPath = "$GOSlimOutDir/input.geneList.gaf";
			&printGeneListGAF($geneGAFOriginalLineHsh_ref, $sgnfcntDEGeneNameHsh_ref, $geneListGAFPath);#->902
			my ($GOColorRamiGOHsh_ref, $geneBasedOneLinerGOSlimResultPath, $GOBasedOneLinerGOSlimResultPath, $geneGOSlimOneToOneResultPath) = &runAndParseMap2Slim($GOSlimOboPathAry_ref, $geneListGAFPath, $regularGOOboPath, $sgnfcntDEGeneNameHsh_ref, $GOSlimOutDir, $DEResultDataHsh_ref, $comparison, $fullDEInfoHeaderAry_ref);#->1616
			my $RamiGOOutDir = "$GOSlimOutDir/svg/";
			my $svgName = $comparison."_".$cutoff."_GOSlim";
			&generateRamiGOTreesvg($GOColorRamiGOHsh_ref, $RamiGOOutDir, $svgName);#->750
			my $svgPath = "$RamiGOOutDir"."$svgName.svg";
			
			$GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'svgPath'} = $svgPath;
			$GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'geneBasedOneLinerGOSlimResultPath'} = $geneBasedOneLinerGOSlimResultPath;
			$GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'GOBasedOneLinerGOSlimResultPath'} = $GOBasedOneLinerGOSlimResultPath;
			$GOSlimResultPathHsh_ref->{$cutoff}{$comparison}{'geneGOSlimOneToOneResultPath'} = $geneGOSlimOneToOneResultPath;
		}
	}
	
	return $GOSlimResultPathHsh_ref, $GOSlimRootDir;
}

exit;
