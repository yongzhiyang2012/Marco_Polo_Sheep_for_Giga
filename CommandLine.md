# Part1 Assembly
### scythe
#### Version: 0.994 BETA
#### Command Line
```
scythe -a adapter_file.fasta -o reads.rmadapter.fq reads.raw.fq
```
### sickle
#### Version: 1.33
#### Command Line
```
sickle pe -f reads.1.rmadapter.fq -r reads.2.rmadapter.fq -t sanger -o reads.1.trim.fq -p reads.2.trim.fq -s reads.single.fq
```
### SOAPec
#### Version: 2.02
#### Command Line
```
KmerFreq_HA -k 23 -t 64 -p panyang -l read.lst -L 100 
Corrector_HA -k 23 -l 3 -r 50 -t 10 panyang.freq.gz read.lst 
```
### platanus 
#### Version: 1.2.4
#### Command Line
```
platanus assemble -o panyang -m 1500 -t 30 -f small.KS1.1.fq small.KS1.2.fq small.KS2.1.fq small.KS2.2.fq small.KS3.1.fq small.KS3.2.fq small.KS4.1.fq small.KS4.2.fq small.KS5.1.fq small.KS5.2.fq  2>01.platanus.pl.contig.sh.log
platanus scaffold -o panyang -c panyang_contig.fa -b panyang_contigBubble.fa -IP1 small.KS1.1.fq small.KS1.2.fq -IP2 small.KS2.1.fq small.KS2.2.fq -IP3 small.KS3.1.fq small.KS3.2.fq -IP4 small.KS4.1.fq small.KS4.2.fq -IP5 small.KS5.1.fq small.KS5.2.fq -OP6 Large-Fragment.KS-1.1.fq Large-Fragment.KS-1.2.fq -OP7 Large-Fragment.KS-2.1.fq Large-Fragment.KS-2.2.fq -OP8 Large-Fragment.KS-3.1.fq Large-Fragment.KS-3.2.fq -OP9 Large-Fragment.KS-4.1.fq Large-Fragment.KS-4.2.fq -OP10 Large-Fragment.KS01-3.1.fq Large-Fragment.KS01-3.2.fq -OP11 Large-Fragment.KS01-4.1.fq Large-Fragment.KS01-4.2.fq -t 30 2>01.platanus.pl.scaffold.sh.log
platanus gap_close -o panyang -c panyang_scaffold.fa -IP1 small.KS1.1.fq small.KS1.2.fq -IP2 small.KS2.1.fq small.KS2.2.fq -IP3 small.KS3.1.fq small.KS3.2.fq -IP4 small.KS4.1.fq small.KS4.2.fq -IP5 small.KS5.1.fq small.KS5.2.fq -OP6 Large-Fragment.KS-1.1.fq Large-Fragment.KS-1.2.fq -OP7 Large-Fragment.KS-2.1.fq Large-Fragment.KS-2.2.fq -OP8 Large-Fragment.KS-3.1.fq Large-Fragment.KS-3.2.fq -OP9 Large-Fragment.KS-4.1.fq Large-Fragment.KS-4.2.fq -OP10 Large-Fragment.KS01-3.1.fq Large-Fragment.KS01-3.2.fq -OP11 Large-Fragment.KS01-4.1.fq Large-Fragment.KS01-4.2.fq -t 30 2>01.platanus.pl.gapclose.sh.log
```
# Part2 Assessment
### BUSCO
#### Version: 2.0.1
#### Command Line
```
BUSCO.py -c 32 -o <out dir> -in panyang.pep -l mammalia_odb9 -m prot
BUSCO.py -c 32 -o <out dir> -sp human -in panyang.fa -l mammalia_odb9 -m geno --long
```
### CEGMA
#### Version: 2.5
#### Command Line
```
cegma -T 30 -g panyang.fa -o <output prefix>
```
### FRCurve
#### Version: 1.3.0
#### Command Line
```
FRC --genome-size 3000000000 --pe-sam lib.500.bam --mp-sam lib.15000.bam --pe-max-insert 500 --mp-max-insert 15000 --out frc_curve
```
# Part3 SNPs and InDels
### BWA
#### Version: 0.7.12
#### Command Line
```
bwa mem -t 12 -R '@RG ID:<sample id>  SM:<sample id>  LB:<sample id>' panyang.fa lib1.1.fq.gz lib1.2.fq.gz | samtools sort -O bam -T ./ -l 3 -o <sample id>.bam -
samtools rmdup <sample id>.bam <sample id>.rmdup.bam
```
### GATK
#### Version: 3.3
#### Command Line
```
java -Xmx10g -jar GenomeAnalysisTK.jar -R panyang.fa -T RealignerTargetCreator -o <sample id>.intervals -I <sample id>.rmdup.bam
java -Xmx10g -jar GenomeAnalysisTK.jar -R panyang.fa -T IndelRealigner -targetIntervals <sample id>.intervals -o <sample id>.realn.bam -I <sample id>.rmdup.bam
java -Xmx10g -jar GenomeAnalysisTK.jar -nct 12 -R panyang.fa -T HaplotypeCaller -I <sample id>.realn.bam -out_mode EMIT_VARIANTS_ONLY -o <sample id>.vcf
```
# Part4 Annotation
### RepeatMasker
#### Version: 1.323
#### Command Line
```
RepeatMasker -nolow -no_is -norna -parallel 1 -species mammal -gff panyang.fa
```
### RepeatProteinMask
#### Version: 1.36
#### Command Line
```
RepeatProteinMask -engine ncbi -noLowSimple -pvalue 0.0001 panyang.fa
```
### RepeatModeler 
#### Version: 1.0.8
#### Command Line
```
BuildDatabase -name panyang panyang.fa
RepeatModeler -pa 30 -database panyang
RepeatMasker -lib RM*/consensi.fa.classified -pa 30 panyang.fa
```
### TRF
#### Version: 407b
#### Command Line
```
trf panyang.fa 2 7 7 80 10 50 2000 -d -h
```
### BLAST
#### Version: 2.28
#### Command Line
```
blastall -p tblastn -d <db> -i <query> -e 1E-5 -o protein2query.out -a 12
```
### BLAST2GENE
#### Version: 17
#### Command Line
```
blast.parse.pl protein2query.out > protein2query.bp
blast2gene.pl protein2query.bp > protein2query.bl2g
```
### GeneWise
#### Version: 2.41
#### Command Line
```
genewise -u <start> -v <end> -<trev|tfor> -gff query.fa <chr name>.fa > genewise.gff
```
### Augustus
#### Version: 2.5.5
#### Command Line
```
augustus --species=human <chr name>.fa > <chr name>.gff
```
### GenScan
#### Version: Not available
#### Command Line
```
genscan HumanIso.smat <chr name>.fa > <chr name>.out
```
### EVM
#### Version: 1.1.1
#### Command Line
```
evidence_modeler.pl --genome <chr name>.fa --weights weights.txt --gene_predictions ab_initio.gff --protein_alignments homolog.gff > evm.out; EVM_to_GFF3.pl evm.out <chr name> > evm.out.gff
```
# Part5 Synteny
### last
#### Version: 761
#### Command Line
```
lastdb -uNEAR -cR11 ref_db ref.fa
lastal -P48 -m100 -E0.05 ref_db panyang.fa | last-split > query2db.maf
maf-swap query2db.maf | last-split > cattle.panyang.sing.maf # same parameters were applied for the other species
```
# Part6 Evolution analysis
### prank
#### Version: 150803
#### Command Line
```
prank -d=input.fa  -o=output.fa -f=fasta -codon
```
### RAxML
#### Version: 7.2.8
#### Command Line
```
raxmlHPC -# 100 -b 12345 -f a -m GTRGAMMAI -s align.phy    # the generated alignments were applied for bootstrap
```
### orthoMCL
#### Version: 2.0.9
#### Command Line
```
orthomclInstallSchema orthomcl.config.template
orthomclAdjustFasta compliantFasta/opa Opossum.pep 1
orthomclAdjustFasta compliantFasta/dog dog.pep 1
orthomclAdjustFasta compliantFasta/hor horse.pep 1
orthomclAdjustFasta compliantFasta/hum human.pep 1
orthomclAdjustFasta compliantFasta/she sheep.pep 1
orthomclAdjustFasta compliantFasta/oam panyang.pep 1
orthomclAdjustFasta compliantFasta/cow cattle.pep 1
orthomclAdjustFasta compliantFasta/pig pig.pep 1
orthomclAdjustFasta compliantFasta/goa goat.pep 1
orthomclFilterFasta compliantFasta/ 10 20
makeblastdb -in goodProteins.fasta -dbtype prot
blastp -db goodProteins.fasta -query goodProteins.fasta -out all-all.blastp.out -evalue 1e-5 -outfmt 6 -num_threads 24
orthomclBlastParser all-all.blastp.out compliantFasta > similarSequences.txt
perl -p -i -e 's/0\t0/1\t-181/' similarSequences.txt
orthomclLoadBlast orthomcl.config.template similarSequences.txt
orthomclPairs orthomcl.config.template orthomcl_pairs.log cleanup=no
orthomclDumpPairsFiles orthomcl.config.template
mcl mclInput --abc -I 1.5 -o mclOutput
orthomclMclToGroups cluster 1 < mclOutput > groups.txt
```
### PAML
#### Version: 4.8
#### Command Line
```
codeml codeml.ctl
```
We used the Codeml program from the PAML package with a branch-site model (runmode = -2, model = 2, NSsites = 2) to detect positively selected genes in focal lineages. A likelihood ratio test was constructed to compare a model that allows sites to be under positive selection on the foreground branch with the null model in which sites may evolve neutrally and under purifying selection. The p-values were computed based on the Chi-square statistic adjusted by the FDR method and genes with adjusted p-value < 0.05 were treated as candidates for positive selection.



