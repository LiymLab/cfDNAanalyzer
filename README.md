# cfDNAanalyzer
cfDNAanalyzer (<ins>c</ins>ell-<ins>f</ins>ree <ins>DNA</ins> sequencing data <ins>analyzer</ins>) is a toolkit for cfDNA whole-genome sequencing data analysis which includes two main parts: 
(1) the extraction of genomic and fragmentatomic features at whole-genome or genomic-region levels;
(2) the processing of extracted features and the building of machine learning models for disease detection and classification.<br>
<summary><h2>Table of Contents</h2></summary>
<li>
  <a href="#Description">Description</a>
  <ul>
    <li><a href="#Environment-requirement-and-installation">Environment requirement and installation</a></li>
    <li><a href="#Supported-features">Supported features</a></li>
    <li><a href="#Supported-feature-processing-methods-and-machine-learning-models">Supported feature processing methods and machine learning models</a></li>
    <li><a href="#Usage">Usage</a></li>
    <li><a href="#Run-the-usage-example">Run the usage example</a></li>
  </ul>
</li>
<li>
  <a href="#Output-files">Output files</a>
  <ul>
    <li><a href="#Features">Features</a></li>
    <li><a href="#Disease-detection-and-classification">Disease detection and classification</a></li>
  </ul>
</li>
<li>
  <a href="#Versions-of-packages-in-our-environment">Versions of packages in our environment</a>
  <ul>
    <li><a href="#Python">Python</a></li>
    <li><a href="#R">R</a></li>
  </ul>
</li>
<li>
  <a href="#Contact">Contact</a>
</li>

## Description
### Environment requirement and installation
Please ensure [<ins>samtools (v1.3.1)</ins>](https://github.com/samtools/samtools), [<ins>bedtools (v2.29.2)</ins>](https://bedtools.readthedocs.io/en/latest/index.html), and [<ins>deeptools (3.5.1)</ins>](https://github.com/deeptools/deepTools) are in your environment. Then, you can install the toolkit following the steps below:
```
git clone https://github.com/LiymLab/cfDNAanalyzer.git
cd cfDNAanalyzer/
conda create -n cfDNAanalyzer --clone ./envs/cfDNAanalyzer
conda activate cfDNAanalyzer
Rscript install_R_packages.R
```
### Supported features
### 1.  Features for whole genome:
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1  <ins>C</ins>opy <ins>N</ins>umber <ins>A</ins>lterations (CNA) ([<ins>Adalsteinsson *et al, Nat Commun*, 2017</ins>](https://www.nature.com/articles/s41467-017-00965-y))<br>
* Copy number alterations comprise deletions or amplifications of a particular region of the genome, with a size as low as a few kilobases up to entire chromosomes. 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2  <ins>E</ins>nd <ins>M</ins>otif frequency and diversity (EM) ([<ins>Zhou et al, 2023</ins>](https://www.pnas.org/doi/10.1073/pnas.2220982120))
* End motifs are the terminal n-nucleotide sequence (4-mer end motif in this toolkit) at each 5′ fragment end of cfDNA molecules. End motif frequency refers to the frequency of all 256 4-mer end motifs.<br>
* End motif diversity is the normalized Shannon entropy of the categorical distribution of all possible 4-mer end-motifs of all cfDNA fragments.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3  <ins>F</ins>ragmentation <ins>P</ins>rofile (FP) ([<ins>Cristiano *et al, Nature*, 2019</ins>](https://doi.org/10.1038/s41586-019-1272-6))
* Fragmentation profile describes fragmentation patterns of cfDNA across the genome, which is the fraction of short cfDNA fragments (100–150 bp) to long cfDNA fragments (151–220 bp) for each 5Mb window across the genome.

### 2.  Features for specific regions: 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.1  <ins>N</ins>ucleosome <ins>O</ins>ccupancy and <ins>F</ins>uzziness (NOF) ([<ins>Li *et al, Genome Med*, 2024</ins>](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01280-6))
* Nucleosome occupancy reflects the frequency with which nucleosomes occupy a given genomic region in a cell population. Nucleosome occupancy for a specific region is calculated as the average occupancy values of all based in this region.<br>
* Nucleosome fuzziness is defined as the deviation of nucleosome positions within a region in a cell population and could reflect cell heterogeneity at the chromatin level. Nucleosome fuzziness for a specific region is calculated as the average fuzziness of all the nucleosomes whose center is located within the region.<br>
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.2 <ins>N</ins>ucleosome <ins>P</ins>rofile (NP) ([<ins>Doebley *et al, Nat Commun*, 2022</ins>](https://www.nature.com/articles/s41467-022-35076-w))
* Nucleosome profile is a composite coverage profile computed as the mean of the GC-corrected cfDNA fragment midpoint coverage across a set of sites (Binding sites sets of 377 transcription factors as the dafalut). For each set of sites, three features are identified from the composite coverage profile:
   * (1) The average coverage value from ± 30 bp of the central site of each set (central coverage).
   * (2) The average coverage value from ± 1000 bp of the central site of each set (mean coverage).
   * (3) The overall nucleosome peak amplitude is calculated by using a Fast Fourier Transform on the window ± 960 bp from the site and taking the amplitude of the 10th frequency term (amplitude).
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.3  <ins>W</ins>indowed <ins>P</ins>rotection <ins>S</ins>core (WPS) ([<ins>Snyder *et al, Cell*, 2016</ins>](https://doi.org/10.1016/j.cell.2015.11.050))
* Windowed protection score is the number of DNA fragments completely spanning a window centered at a given genomic coordinate minus the number of fragments with an endpoint within that same window. WPS can be calculated using long fragments (120–180 bp; 120 bp window) or short fragments (35–80 bp; 16 bp window). The WPS for a specific region is defined as the average WPS of all bases in this region. High WPS values indicate increased protection of DNA from digestion while low values indicate that DNA is unprotected.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.4  <ins>O</ins>rientation-aware <ins>C</ins>fDNA <ins>F</ins>ragmentation (OCF) ([<ins>Sun *et al, Genome Res.*, 2019</ins>](https://genome.cshlp.org/content/29/3/418.long))
* Orientation-aware cfDNA fragmentation is the differences of read densities of the upstream and downstream fragment ends in specific genomic regions.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.5  <ins>E</ins>nd <ins>M</ins>otif frequency and diversity for <ins>R</ins>egions (EMR)([<ins>Zhou et al, 2023</ins>](https://www.pnas.org/doi/10.1073/pnas.2220982120))
* We introduced end motif frequency and diversity for regions, which is defined as the frequency and diversity of all 256 4-mer end motifs for each region.<br>
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.6  <ins>F</ins>ragmentation <ins>P</ins>rofile for <ins>R</ins>egions (FPR) ([<ins>Cristiano *et al,Nature*, 2019</ins>](https://doi.org/10.1038/s41586-019-1272-6))
* We introduced fragmentation profile for regions, which is defined as the fraction of short cfDNA fragments (100–150 bp) to long cfDNA fragments (151–220 bp) for each region.

### 3.  Features for transcription start sites (TSSs): 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.1 <ins>P</ins>romoter <ins>F</ins>ragmentation <ins>E</ins>ntropy (PFE) ([<ins>Esfahani *et al, Nat Biotechnol*, 2022</ins>](https://doi.org/10.1038/s41587-022-01222-4))<br>
* Promoter fragmentation entropy quantifies the diversity of cfDNA fragment lengths around the TSSs of genes. It is calculated by a modified Shannon index for lengths of cfDNA fragment where both ends fell within ±1 kb of the TSS. Then this cfDNA entropy measure is adjusted using a Dirichlet-multinomial model for normalization.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2 <ins>TSS</ins> <ins>C</ins>overage (TSSC) ([<ins>Ulz *et al, Nat Genet*, 2016</ins>](https://www.nature.com/articles/ng.3648))
* TSS coverage refers to the cfDNA sequencing read coverage in the region surrounding TSSs.

### Supported feature processing methods and machine learning models

### Usage
```ruby
bash cfDNAanalyzer.sh -I <InputFile> -o <OutputDirectory> -F <Features> [Options]
``` 
#### Options: 
```
-- General options
  -I  FILE      A text file containing all input BAM files with one BAM file per line. 
                BAM files generated using both Bowtie2 and BWA are accepted.  
  -o  DIR       Output directory for all the results. Default: [./]
  -F  STR       Features to extract, including CNA, NOF, WPS, EM, EMR, FP, FPR, NP, OCF, PFE, and TSSC. 
                Features should be set as a string separated by comma, e.g., CNA,NOF. 
                Default: All available features will be extracted.
                Note: The following features are specifically designed for paired-end sequencing data: FP, FPR, NP, PFE, and OCF.
  -g  STR       Genome version of input BAM files (hg19/hg38). Default: [hg38] 
  -b  FILE      A BED3 file specifying the regions to extract features.
                The file should contain three TAB-delimited columns: chromosome start end.              

-- Options specific for Copy Number Alterations (CNA)
  -B  INT       Bin size in kilobases (10, 50, 500, or 1000). Default: [1000]
  --CNA  STR    Additional parameter setting for software ichorCNA. 
                The full parameter list is available by running xxx. [optional]

-- Options specific for Nucleosome Occupancy and Fuzziness (NOF)
  --NOF  STR    Additional parameter setting for software DANPOS2. 
                The full parameter list is available by running xxx. [optional]

-- Options specific for Windowed Protection Score (WPS)
  -x  INT       Min fragment length used for long fragments WPS calculation. Default: [120]
  -X  INT       Max fragment length used for long fragments WPS calculation. Default: [180]
  -w  INT       Window size for long fragments WPS calculation. Default: [120]
  -m  INT       Min fragment length used for short fragments WPS calculation. Default: [35]
  -M  INT       Max fragment length used for short fragments WPS calculation. Default: [80]
  -W  INT       Window size for short fragments WPS calculation. Default: [16]

-- Options specific for End Motif frequency and diversity (EM)
  -f  FILE      Reference genome in FASTA format. 
                For example, hg38 reference genome can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz.

-- Options specific for Nucleosome Profile (NP)
  -l  DIR       Directory containing a list of files with each file for a set of sites.
                The file must have at least two columns with the following column names: "Chrom" for the chromosome name and "position" for the site positions. If not provided, the 377 TF binding site lists from the referenced Nucleosome Profile paper will be used (xxx). 

-- Options for Promoter Fragmentation Entropy (PFE)
  --PFE  STR    Addtional parameter setting for PFE analysis. The full parameter list is available by running xxx.[optional]

-- Options for TSS Coverage (TSSC)
  -u  INT                    Number of base pairs upstream of TSS used for calculating TSSC. Default: [1000]
  -d  INT                    Number of base pairs downstream of TSS used for calculating TSSC. Default: [1000]
  -S  FILE                   A BED6 file specifying the coordinates of TSSs used for calculating TSSC.  
  --bamCoverage  STR         Additional parameter seting for the "bamCoverage" command of deeptools. [optional] 
  --multiBigwigSummary  STR  Additional parameter seting for the "multiBigwigSummary" command of deeptools. [optional]
```
                        
### Run the usage example
You can directly run cfDNAanalyzer by ```cfDNAanalyzer.sh``` provided in the ```cfDNAanalyzer/``` directory with the example files in xxx.
```
bash cfDNAanalyzer.sh -I ./input/bam_input.txt -o ./output/ -F CNV,NOF,TSS,WPS,EM,FP,NP,PFE,OCF -g hg19 -f <Reference.fa> -b ./End_motif_frequency/tss_2k_regions.bed > ./cfDNAanalyzer.log
``` 

## Output files
### Features
#### Copy Number Alterations (CNA)
```CNV.txt``` is xxx with four columns specifying the chromosome, start coordinate, end coordinate, and the estimated copy number for each bin. 
```r
chr	start	end	sample.copy.number
1	1000001	2000000	2
1	2000001	3000000	2
1	4000001	5000000	2
```
```CNV.wig``` stores xxx.<br>
```r
fixedStep chrom=chr1 start=1 step=1000000 span=1000000
1792
4523
5506
```

#### Nucleosome Occupancy and Fuzziness (NOF)
```occupancy.bed``` stores the nucleosome occupancy features with five columns specifying the chromosome, start coordinate, end coordinate, region ID, and occupany value for each region.
```r
chr	start	end	ID	occupancy
chr1	68091	70091	region_1	6.616
chr1	138379	140379	region_2	19.206
chr1	366640	368640	region_3	2.88
```
&nbsp;<br>
```meanfuziness.tsv``` stores the nucleosome fuzziness features with four columns specifying the chromosome, start coordinate, end coordinate, and the average fuzziness value for each region.
```r
chr  start  end  meanfuziness
chr1 68091 70091 33.4683
chr1 138379 140379 48.2868
chr1 366640 368640 37.744
```
&nbsp;<br>
```pooled/<sampleID>.smooth.wig``` stores the occupancy values across the entire genome for each processed BAM file.<br>
```r
fixedStep chrom=chr1 start=1  step=10 span=10
0.0
3.0
5.2
```
&nbsp;<br>
```pooled/<sampleID>.smooth.positions.xls``` store the fuzziness values of all identified nucleosomes in each processed BAM file. This file has six columns specifying the chromosome, start coordinate, end coordinate, summit position, occupancy value at the summit position, and fuzziness score for each nucleosome.
```r
chr     start   end     smt_pos smt_value       fuzziness_score
chr1    10011   10151   10081   24.0    56.2199218395978
chr1    10161   10301   10231   56.0    52.17928320347469
chr1    10301   10441   10371   48.0    52.888443917362515
```

#### Windowed Protection Score (WPS)
```WPS.txt``` stores the WPS values with five columns specifying the chromosome, start coordinate, end coordinate, WPS for long fragments, and WPS for short fragments for each region.
```r
chr    start    end    long_WPS    short_WPS
chr1	68091	70091	-0.781609	0
chr1	138379	140379	-1.94753	0.0164918
chr1	366640	368640	-0.338331	0
```

#### End Motif frequency and diversity for regions (EMR)
Our toolkit outputs two types of EMR features at different levels:
(1) Motif frequency and diversity for all input regions aggregated together.
```all_motifs_frequency.txt```<br>
```r
Motif   Frequency
CTAT    0.0030392622975323673
ATAG    0.0030392622975323673
TATT    0.0059179091631653994
```
```all_motifs_mds.txt```<br>
```r
MDS     0.9747882739505999
```
&nbsp;<br>
(2) Motif frequency and diversity for each region separately. "index" refers to the line number of the region. 
```region_<index>_motif_frequency_and_mds.txt```<br>
```r
Motif   Frequency
CTAT    0.006259389083625438
ATAG    0.006259389083625438
TATT    0.009764646970455683
```
#### End Motif frequency and diversity of whole genome (EM)

#### Fragmentation Profile for regions (FPR)
```Fragmentation_Profile.txt``` has six columns specifying the chromosome, start coordinate, end coordinate, number of short fragments, number of long fragments, and ratio of short to long fragments for each input region.
```r
seqnames	start	end	short	long	ratio
chr1	738137	740137	1	3	0.333333333333333
chr1	817043	819043	5	5	1
chr1	865445	867445	1	3	0.333333333333333
```

#### Fragmentation Profile of whole genome (FP)

#### Nucleosome Profile (NP)
```NucleosomeProfile.txt``` has four columns specifying the file name of sites set, the mean coverage, central coverage, and nucleosome peak amplitude for each sites set in the input directory. 
```r
site_name	mean_coverage	central_coverage	amplitude
IKZF1.hg38.10000.txt	0.99428	1.01050	0.45213
TCF7L1.hg38.10000.txt	1.00737	1.03907	0.20425
NKX3-1.hg38.10000.txt	0.99402	1.00351	0.27422
```
```plots/<site_list>.pdf``` stores the cfDNA fragment coverage profile for each sites set in the input directory. <br>

#### Orientation-aware CfDNA Fragmentation (OCF)
```OCF.txt```<br>
The "TSS_ID" column is the id of each TSS 2k region.<br>
The "OCF_Ratio" column is orientation-aware cfDNA fragmentation ratio for each TSS 2k region.<br>
```r
TSS_ID	OCF_Ratio
A1BG_1	-0.142857142857143
A1CF_1	-0.2
A2M_1	0.333333333333333
```

#### Promoter Fragmentation Entropy (PFE)
```PFE.txt``` has two columns specifying the TSS ID (The fourth column in the input file?) and PFE value for the 2kb surrounding TSS.
```r
TSS_ID	PFE
A1BG_1	0.341049657283462
A1CF_1	0.34067802945418
A2M_1	0.340071343470199
```

### TSS Coverage (TSSC)
```average_coverage.txt``` has four columns specifying the chromosome, start coordinate, end coordinate, and cfDNA fragment coverage value for each TSS surrounding regions. 
```r
chr	start	end	coverage
chr1	68091	70091	1.8115
chr1	138379	140379	5.7315
chr1	366640	368640	0.79
```

### Disease detection and classification


## Versions of packages in our environment
### Python:
```r
python                         3.7.16

alembic                        1.8.1
anyio                          3.5.0
appdirs                        1.4.4
argon2-cffi                    20.1.0
async-generator                1.10
attrs                          22.1.0
Babel                          2.11.0
backcall                       0.2.0
backports.zoneinfo             0.2.1
beautifulsoup4                 4.11.1
bio                            1.6.2
biopython                      1.81
biothings-client               0.3.1
bleach                         4.1.0
blinker                        1.4
Boruta                         0.3
Bottleneck                     1.3.5
brotlipy                       0.7.0
bwa                            1.1.1
bx                             0.3.0
bx-python                      0.10.0
cachetools                     5.3.3
certifi                        2023.7.22
certipy                        0.1.3
cffi                           1.15.1
charset-normalizer             2.0.4
conda-package-handling         2.0.2
conda_package_streaming        0.7.0
ConfigArgParse                 1.7
cryptography                   39.0.1
cssselect                      1.2.0
cssutils                       2.7.1
cycler                         0.11.0
datrie                         0.8.2
debugpy                        1.5.1
decorator                      5.1.1
deepTools                      3.5.1
deeptoolsintervals             0.1.9
defusedxml                     0.7.1
docutils                       0.20.1
entrypoints                    0.4
exceptiongroup                 1.2.0
fastjsonschema                 2.16.2
fonttools                      4.25.0
gitdb                          4.0.11
GitPython                      3.1.43
gprofiler-official             1.0.0
greenlet                       2.0.1
h11                            0.14.0
httpcore                       0.17.3
httpx                          0.24.1
idna                           3.4
importlib-metadata             4.11.3
importlib-resources            5.2.0
iniconfig                      2.0.0
ipykernel                      6.15.2
ipython                        7.31.1
ipython-genutils               0.2.0
ipywidgets                     7.6.5
jedi                           0.18.1
Jinja2                         3.1.2
joblib                         1.3.0
json5                          0.9.6
jsonschema                     4.17.3
jupyter                        1.0.0
jupyter_client                 7.4.9
jupyter-console                6.4.4
jupyter_core                   4.11.2
jupyter-server                 1.23.4
jupyter-telemetry              0.1.0
jupyterhub                     3.1.1
jupyterhub-dummyauthenticator  0.3.1
jupyterlab                     3.5.3
jupyterlab-language-pack-zh-CN 4.0.post0
jupyterlab-pygments            0.1.2
jupyterlab_server              2.19.0
jupyterlab-widgets             1.0.0
kiwisolver                     1.4.4
lightgbm                       4.3.0
lxml                           5.2.1
Mako                           1.2.3
MarkupSafe                     2.1.1
matplotlib                     3.4.1
matplotlib-inline              0.1.6
mHapTk                         1.0
mistune                        0.8.4
mkl-fft                        1.3.1
mkl-random                     1.2.2
mkl-service                    2.4.0
mlxtend                        0.23.1
modules                        1.0.0
munkres                        1.1.4
mygene                         3.2.2
nbclassic                      0.5.2
nbclient                       0.5.13
nbconvert                      6.4.4
nbformat                       5.7.0
nest-asyncio                   1.5.6
nose                           1.3.7
notebook                       6.5.2
notebook_shim                  0.2.2
numexpr                        2.8.4
numpy                          1.21.6
oauthlib                       3.2.0
packaging                      22.0
pamela                         1.0.0
pandas                         1.3.2
pandocfilters                  1.5.0
parso                          0.8.3
pexpect                        4.8.0
pickleshare                    0.7.5
Pillow                         9.4.0
pip                            22.3.1
pkgutil_resolve_name           1.3.10
platformdirs                   4.0.0
plotly                         5.9.0
pluggy                         1.0.0
ply                            3.11
pooch                          1.8.1
premailer                      3.10.0
prometheus-client              0.14.1
prompt-toolkit                 3.0.36
psutil                         5.9.0
ptyprocess                     0.7.0
py2bit                         0.3.0
pybedtools                     0.8.0
pyBigWig                       0.3.17
pycosat                        0.6.4
pycparser                      2.21
Pygments                       2.11.2
PyJWT                          2.4.0
pyOpenSSL                      23.0.0
pyparsing                      3.0.9
PyQt5-sip                      12.11.0
pyrsistent                     0.18.0
pysam                          0.19.0
PySocks                        1.7.1
pytest                         7.4.4
python-dateutil                2.8.2
python-json-logger             2.0.1
python-telegram-bot            20.3
pytz                           2022.7
PyYAML                         3.12
pyzmq                          23.2.0
qtconsole                      5.4.0
QtPy                           2.2.0
ratelimiter                    1.2.0.post0
requests                       2.28.1
rpy2                           3.3.3
ruamel.yaml                    0.17.21
ruamel.yaml.clib               0.2.6
scikit-learn                   1.0.2
scipy                          1.7.1
Send2Trash                     1.8.0
setuptools                     65.6.3
sip                            6.6.2
six                            1.16.0
sklearn                        0.0.post5
smmap                          5.0.1
snakemake                      5.19.2
sniffio                        1.2.0
soupsieve                      2.3.2.post1
SQLAlchemy                     1.4.39
Tcl                            0.2
tenacity                       8.0.1
terminado                      0.17.1
testpath                       0.6.0
threadpoolctl                  3.1.0
toml                           0.10.2
tomli                          2.0.1
toolz                          0.12.0
toposort                       1.10
torch                          1.13.1
torchaudio                     0.13.1
torchvision                    0.14.1
tornado                        6.2
tqdm                           4.64.1
traitlets                      5.7.1
typing_extensions              4.7.1
tzlocal                        5.1
urllib3                        1.26.14
wcwidth                        0.2.5
webencodings                   0.5.1
websocket-client               0.58.0
wheel                          0.38.4
widgetsnbextension             3.5.2
wrapt                          1.16.0
xgboost                        1.6.2
yagmail                        0.15.293
zipp                           3.11.0
zstandard                      0.19.0
```

### R:
```r
R                              4.3.0

DescTools                      0.99.40
zoo                            1.8.12
plyr                           1.8.9
reshape2                       1.4.4
data.table                     1.15.2
MASS                           7.3-60.0.1
e1071                          1.7-14
gtools                         3.9.5
matrixStats                    1.2.0
optparse                       1.7.4
httr                           1.4.7
tidyverse                      2.0.0
RCurl                          1.98-1.14
BiocManager                    3.18.1
HMMcopy                        1.44.0
GenomeInfoDb                   1.38.8
GenomicRanges                  1.54.1
Rsamtools                      2.18.0
GenomicAlignments              1.38.2
biovizBase                     1.50.0
BSgenome.Hsapiens.UCSC.hg19    1.4.3
BSgenome.Hsapiens.UCSC.hg38    1.4.5
```
## Contact
Yumei Li: ymli12@suda.edu.cn
