# Sequence diversity calculation -- Guanqiao Feng 2019/05/20

## Step 1. Extract sites from bam files
which meet these criteria: (1) mapping quality >= 20, (2) read quality >= 20; (3) 10-150 coverage; (4) >=75% genotype

For each \*dedup.bam file:
```sh
$ samtools depth -q 20 -Q 20 /home/OlsonLab_Data/seq_cap/salixSDR/goodingii/alignments/SN503M.dedup.bam > SN503M.depth &
# extract depth info for sites with (1) mapping quality >=20 and (2) read quality >= 20

$ awk -F "\t" '{ if (($3>=10) && ($3<=150)) { print } }' SN503M.depth | cut -f 1,2 > SN503M_10TO150.depth &
# extract sites with (3) depth 10-150
```
Generate a folder for male and female each, calculating male and female sequence diversity seperately, using male as an example. Codes are the same for females
```sh
cat *depth | cut -f 1,2 | sort | uniq -c > SN_M_10TO150_count.depth &
# merge sites which matches criteria (1)(2)(3) for all males

sed 's/ \+/\t/g' SN_M_10TO150_count.depth > SN_M_10TO150_count_tab.depth
# reformat the merged file

awk -F "\t" '{ if ($2 >=18) { print $2"\t"$3"\t"$4 } }' SN_M_10TO150_count_tab.depth > SN_M_10TO150_count_tab_025cov.depth
# screen with (4) >=75% genotype. For S. nigra male, total individual number is 24. 75% of 24 is 18. Change the value 18 accordingly if you total number is not 24

cut -f 2,3 SN_M_10TO150_count_tab_025cov.depth | sort -k1,1 -k2n > SN_M_10TO150_count_tab_025cov_pos.depth
# reformat the site file
```

## Step 2. Calculate diversity with vcftools
```sh
vcftools --vcf /home/gfeng/S_nigra_variant/depth/snigra.genotypes.snp.vcf --site-pi --positions SN_M_10TO150_count_tab_025cov_pos.depth --keep males.txt --out SN_M_SNP.txt
# it calculate pi for each site with the given sites identified from step 1 with the given individuals privided by --keep
# the output file contains two types of results 1) non-zero pi for heterozygosity sites 2) zero pi for homozygosity sites which are different from the reference
# we will use the non-zero pi
```

## Step 3. Map diversity values (from step 2) to the sites (from step 1) and calculate diversity
```sh
perl add_diversity.pl
# it maps diversity values to the sites

perl per_site_diversity.pl
# it calculate per site diversity using 1000bp as window

perl extract_SDR.perl
# it extracted sex determination regions

perl per_site_diversity_SDR.pl
# it calculate per site diversity in sex determination regions
```
In the end, copy all *_persite.depth files (3 from male and 3 from female folder) to a folder which we are going to use to generate figures

### Step 4. Figure generation
Generate the figure with R:
```sh
setwd("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/code")

male_genome <- read.csv("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/data/SN_M_10TO150_count_tab_025cov_pos_diversity_persite.depth", head = T, sep = "\t")
male_genome_pi = male_genome$pi_per_site
male_genome_theta = male_genome$theta_per_site

male_SNsdr <- read.csv("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/data/SN_M_10TO150_count_tab_025cov_pos_diversity_SNsdr_persite.depth", head = T, sep = "\t")
male_SNsdr_pi = male_SNsdr$pi_per_site
male_SNsdr_theta = male_SNsdr$theta_per_site

male_SEsdr <- read.csv("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/data/SN_M_10TO150_count_tab_025cov_pos_diversity_SEsdr_persite.depth", head = T, sep = "\t")
male_SEsdr_pi = male_SEsdr$pi_per_site
male_SEsdr_theta = male_SEsdr$theta_per_site

female_genome <- read.csv("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/data/SN_F_10TO150_count_tab_025cov_pos_diversity_persite.depth", head = T, sep = "\t")
female_genome_pi = female_genome$pi_per_site
female_genome_theta = female_genome$theta_per_site

female_SNsdr <- read.csv("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/data/SN_F_10TO150_count_tab_025cov_pos_diversity_SNsdr_persite.depth", head = T, sep = "\t")
female_SNsdr_pi = female_SNsdr$pi_per_site
female_SNsdr_theta = female_SNsdr$theta_per_site

female_SEsdr <- read.csv("/Users/gfeng/Documents/project/LAB/S_nigra/diversity/data/SN_F_10TO150_count_tab_025cov_pos_diversity_SEsdr_persite.depth", head = T, sep = "\t")
female_SEsdr_pi = female_SEsdr$pi_per_site
female_SEsdr_theta = female_SEsdr$theta_per_site



library(vioplot)
pdf(width = 10, height = 8, "/Users/gfeng/Documents/project/LAB/S_nigra/diversity/figure/Pi.pdf")
vioplot(female_genome_pi, male_genome_pi, female_SNsdr_pi, male_SNsdr_pi, female_SEsdr_pi, male_SEsdr_pi, names = c("F_genome", "M_genome", "F_SN_SDR", "M_SN_SDR", "F_SE_SDR", "M_SE_SDR"), col = c("orange", "skyblue", "orange", "skyblue", "orange", "skyblue"))
title("Pi of S. nigra")
dev.off()

pdf(width = 10, height = 8, "/Users/gfeng/Documents/project/LAB/S_nigra/diversity/figure/Theta.pdf")
vioplot(female_genome_theta, male_genome_theta, female_SNsdr_theta, male_SNsdr_theta, female_SEsdr_theta, male_SEsdr_theta, names = c("F_genome", "M_genome", "F_SN_SDR", "M_SN_SDR", "F_SE_SDR", "M_SE_SDR"), col = c("orange", "skyblue", "orange", "skyblue", "orange", "skyblue"))
title("Theta of S. nigra")
dev.off()

pdf(width = 16, height = 8, "/Users/gfeng/Documents/project/LAB/S_nigra/diversity/figure/Pi_Theta.pdf")
vioplot(female_genome_theta, female_genome_pi, male_genome_theta, male_genome_pi, female_SNsdr_theta, female_SNsdr_pi, male_SNsdr_theta, male_SNsdr_pi, female_SEsdr_theta, female_SEsdr_pi, male_SEsdr_theta, male_SEsdr_pi, names = c("genome_theta", "genome_pi", "genome_theta", "genome_pi", "SNSDR_theta", "SNSDR_pi", "SNSDR_theta", "SNSDR_pi", "SESDR_theta", "SESDR_pi", "SESDR_theta", "SESDR_pi"), col = c("orange", "orange","skyblue", "skyblue","orange", "orange","skyblue", "skyblue","orange", "orange","skyblue","skyblue"))
title("Theta and Pi of S. nigra")
dev.off()


t.test(female_SNsdr_theta, male_SNsdr_theta)

t.test(female_SNsdr_pi, male_SNsdr_pi)
```
