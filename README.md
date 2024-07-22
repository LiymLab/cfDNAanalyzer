# cfDNAanalyzer
This is an integrated tool for extracting the genetic or epigenetic features of cfDNA from a bam file. These features can be extracted <br>
&nbsp;1.for the whole genome, like <ins>C</ins>opy <ins>N</ins>umber <ins>V</ins>ariation (CNV);<br>
&nbsp;2.for given regions, like <ins>N</ins>ucleosome <ins>O</ins>ccupancy and <ins>F</ins>uzziness (NOF), <ins>W</ins>indowed <ins>P</ins>rotection <ins>S</ins>core (WPS), <ins>E</ins>nd <ins>M</ins>otif frequency and diversity (EM),<ins>F</ins>ragmentation <ins>P</ins>rofile (FP), <ins>N</ins>ucleosome <ins>P</ins>rofile (NP), and <ins>O</ins>rientation-aware <ins>C</ins>fDNA <ins>F</ins>ragmentation (OCF);<br>
&nbsp;3.for specific Transcription Start Sites (TSS), like <ins>P</ins>romoter <ins>F</ins>ragmentation <ins>E</ins>ntropy (PFE) and <ins>TSS</ins> coverage (TSSC).<br> 
For the extracted features, this tool also offers a customized pipeline for downstream applications, like cancer detection or tumor subtype classification. 

## Description
### Features can be extracted from this tool
#### 1.Features extracted for the whole genome
&nbsp;&nbsp;&nbsp;1.1Copy Number Variation([<ins>Adalsteinsson et al, 2017</ins>](https://www.nature.com/articles/s41467-017-00965-y))
&nbsp;&nbsp;&nbsp;Copy number variation is the number of copy variation in a region.

#### 2.Features extracted for given regions
##### <ins>N</ins>ucleosome <ins>O</ins>ccupancy and <ins>F</ins>uzziness (NOF) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>https://www.researchsquare.com/article/rs-1889126/v1</ins>](https://www.researchsquare.com/article/rs-1889126/v1)
Nucleosome occupancy, which reflects the frequency with which nucleosomes occupy a given DNA region in a cell population.<br>
Nucleosome fuzziness, which is defined as the deviation of nucleosome positions within a region in a cell population and could reflect cell heterogeneity at the chromatin level.

##### <ins>W</ins>indowed <ins>P</ins>rotection <ins>S</ins>core (WPS) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>http://dx.doi.org/10.1016/j.cell.2015.11.050</ins>](http://dx.doi.org/10.1016/j.cell.2015.11.050)
A per-base WPS is calculated by subtracting the number of fragment endpoints within a window from the number of fragments completely spanning the window. High WPS values indicate increased protection of DNA from digestion; low values indicate that DNA is unprotected.

##### <ins>E</ins>nd <ins>M</ins>otif frequency and diversity (EM) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>https://www.pnas.org/doi/10.1073/pnas.2220982120</ins>](https://www.pnas.org/doi/10.1073/pnas.2220982120)
End motifs were determined from the terminal 4-nucleotide sequence, i.e., 4-mer end motif, at each 5′ fragment end of cfDNA molecules. End motif frequency of each of the motifs (i.e., a total of 256 motifs) was determined from the total number of fragment ends.<br>
End motif diversity is the normalized Shannon entropy of the categorical distribution of all possible end-motif k-mers.

##### <ins>F</ins>ragmentation <ins>P</ins>rofile (FP) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>https://doi.org/10.1038/s41586-019-1272-6</ins>](https://doi.org/10.1038/s41586-019-1272-6)
Fragmentation profile is number of short reads(100–150 bp), number of long reads(151–220 bp) and short/long ratio in a region.

##### <ins>N</ins>ucleosome <ins>P</ins>rofile (NP) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>https://www.nature.com/articles/s41467-022-35076-w</ins>](https://www.nature.com/articles/s41467-022-35076-w)
Nucleosome profile is the distribution of nucleosomes in extracellular DNA in a site list, and we extracted 3 features from each coverage profile. "central coverage" is the coverage value from ± 30 bp of central site. "mean coverage" is the coverage value from ± 1000 bp of central site. The amplitude of the nucleosome peaks surrounding the central site is calculated by using a Fast Fourier Transform on the window ±960 bp from the central site.

##### <ins>P</ins>romoter <ins>F</ins>ragmentation <ins>E</ins>ntropy (PFE) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>https://doi.org/10.1038/s41587-022-01222-4</ins>](https://doi.org/10.1038/s41587-022-01222-4)
Promoter fragmentation entropy is the Shannon entropy of fragments around the ± 1000 bp of transcription start sites.

##### <ins>O</ins>rientation-aware <ins>C</ins>fDNA <ins>F</ins>ragmentation (OCF) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>http://www.genome.org/cgi/doi/10.1101/gr.242719.118</ins>](http://www.genome.org/cgi/doi/10.1101/gr.242719.118)
Orientation-aware cfDNA fragmentation is the differential phasing of upstream (U) and downstream (D) fragment ends in tissue-specific open chromatin regions.

#### 3.Features extracted for specific Transcription Start Sites
##### <ins>TSS</ins> coverage (TSS) &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; [<ins>https://www.nature.com/articles/ng.3648</ins>](https://www.nature.com/articles/ng.3648)
TSS coverage is the coverage around transcription start sites.
### Environment and installation
Before using this tool, we advise using our packaged conda environment and running the script ```install_R_packages.R``` to avoid package version conflicts and the following errors.<br> 
First, navigate to the directory ```cfDNAanalyzer/``` and run the following commands:
```
conda create -n cfDNAanalyzer --clone ./envs/cfDNAanalyzer
conda activate cfDNAanalyzer
Rscript install_R_packages.R
```

### Tools needed for cfDNAanalyzer
```
bedtools                       2.29.2
conda                          23.1.0
deeptools                      3.5.1
samtools                       1.3.1
```

### Input and output
This tool only accepts the bam file as the input file type. You need to provide a text file containing the path of your bam files in each line. We will create a folder for each bam file in the output directory. Every bam file's single feature will be output to a folder.


### Usage：
```
   bash cfDNAanalyzer.sh -I <InputTxt> -o <OutputDirectory> -F <Features> -f <ReferenceType> -r <Reference.fa> -s <SequencingMethod> -b <RegionBed> [Options]
```
   
### Options: 
```
-- Options for all the features
  -I  File                   Input bam file list. Default: [ ]
  -o  PATH                   Output directory for all results. Default: [./]
  -F                         CfDNA features, including CNV, NOF, TSS, WPS, EM, FP, NP, PFE and OCF; if not provided then will extract all the features. Default: [NULL]
  -f                         Type of reference fasta file, including hg19 and hg38. Default: [hg38]
  -r  File                   Reference fasta file. For example, if you want to use hg38 reference, you can download it from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz. Default: [NULL] 
  -s                         Type of sequencing method, including pair and single. Default: [pair]
  -b  File                   Region bed file. This bed file must be a TAB-delimited bed3 file without any header. if not provided then will use +/- 1kb region from TSS sites. Default: [NULL]


-- Options for copy number variation
  -c  INT                    Total clonal CN states. Default: [7]
  -B  INT                    Number of kb for bin size, including 10, 50, 500, 1000. Default: [1000]
  --addCNV                   Addtional parameters for for copy number variation. Default: [NULL] 

-- Options for nucleosome occupancy and fuzziness
  --addNOF                   Addtional parameters for nucleosome occupancy and fuzziness. Default: [NULL] 

-- Options for windowed protection score
  -x  INT                    Minimum length of reads identified as long reads. Default: [120]
  -X  INT                    Maximum length of reads identified as long reads. Default: [180]
  -w  INT                    Size of window to extract WPS for long reads. Default: [120]
  -m  INT                    Minimum length of reads identified as short reads. Default: [35]
  -M  INT                    Maximum length of reads identified as short reads. Default: [80]
  -W  INT                    Size of window to extract WPS for short reads. Default: [16]
  
-- Options for nucleosome profile
  -l  PATH                   Path of site list. This path must only contain site list for nucleosome profile analysis.
                             And each site file must contain two columns. The chromosome column must have "Chrom"" as the header and the site column must have "position"" as the header; 
                             if not provided then will use 377 site lists provided by the paper of nucleosome profile. Default: [NULL]

-- Options for promoter fragmentation entropy or orientation-aware cfDNA fragmentation
  -T  File                   TSS information files. The file must have five columns: (1) chromosome, (2) 1-base TSS coordinate, 
                             (3) Hugo symbol of the gene corresponding to the TSS, (4) Category 
                             (*negativeControl* must exist, other categories are not important), (5) gene 
                             transcript strand (+1/-1) and (6) a column for TSS ID; 
                             if not provided then will use TSS information file provided by the paper of promoter fragmentation entropy. Default: [NULL]
  -U  INT                    Number of bp upstreams the TSS sites. Default: [1000]
  -D  INT                    Number of bp downstreams the TSS sites. Default: [1000]  
  -p  File                   Either the selector bed file or a no header bed file around the TSS sites in the tssinfo file; 
                             if not provided then will use TSS 2k file provided by the paper of promoter fragmentation entropy. Default: [NULL] (Parameter -U/-D and -p can not set together).
  --addEpicSeq               Addtional parameters for Epic Seq. Default: [NULL] 

-- Options for TSS coverage 
  -u  INT                    Number of bp upstreams the TSS sites. Default: [1000]
  -d  INT                    Number of bp downstreams the TSS sites. Default: [1000]
  -S  File                   Region bed file around transcription start sites. This bed file must be a TAB-delimited bed3 file without any header. if not provided then will use +/- 1kb region from TSS sites. Default: [NULL]
                             (Parameter -u/-d and -S can not be set together)
  --addbamCoverage           Addtional parameters for bamCoverage command of deeptools. Default: [NULL] 
  --addmultiBigwigSummary    Addtional parameters for multiBigwigSummary command of deeptools. Default: [NULL]
```
                        
### Run cfDNAanalyzer
The easiest way to manually run cfDNAanalyzer is to use ```cfDNAanalyzer.sh``` provided in the ```cfDNAanalyzer/``` directory. Here is an example of how to launch the shell script from the command line:
```
bash cfDNAanalyzer.sh -I ./input/bam_input.txt -o ./output/ -F CNV,NOF,TSS,WPS,EM,FP,NP,PFE,OCF -f hg19 -r <Reference.fa> -s pair -b ./End_motif_frequency/tss_2k_regions.bed > ./cfDNAanalyzer.log
``` 


### Note: 
The following features are specific for pair-end bam files:<br>
(1) fragmentation profile<br>
(2) nucleosome profile<br>
(3) promoter fragmentation entropy<br>
(4) orientation-aware cfDNA fragmentation<br>


## Output file for every feature

### Feature copy number variation
```CNV.txt```<br>
"chr" column is the chromosome where each bin is located.<br>
"start" column is the starting site of a bin.<br>
"end" column is the ending site of a bin.<br>
"sample.copy.number" column is the estimated copy number in each bin.<br>
![CNV.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721310992220.png)<br>
<br>
```CNV.wig```<br>
Wig file extracted from input bam file.<br>
![CNV.wig](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311015152.png)<br>
<br>


### Feature nucleosome occupancy and fuzziness
```meanfuziness.tsv```<br>
"chr" column is the chromosome where each region is located.<br>
"start" column is the starting site of a region.<br>
"end" column is the ending site of a region.<br>
"meanfuziness" column is the mean fuziness value for each region.<br>
![meanfuziness.tsv](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311045217.png)<br>
<br>
```occupancy.bed```<br>
"chr" column is the chromosome where each region is located.<br>
"start" column is the starting site of a region.<br>
"end" column is the ending site of a region.<br>
"ID" column is the specific id for each region.<br>
"occupancy" column is the occupancy value for each region.<br>
![occupancy.bed](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311077440.png)<br>
<br>
```pooled/<sampleID>.smooth.positions.xls```<br>
"chr" column is the chromosome where each region is located.<br>
"start" column is the starting site of a region.<br>
"end" column is the ending site of a region.<br>
"smt_pos" column is the occupancy summit point for each region.<br>
"smt_value" column is the occupancy value at smt_pos.<br>
"fuzziness_score" column is the fuzziness score for each region.<br>
![pooled/<sampleID>.smooth.positions.xls](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311120508.png)<br>
<br>
```pooled/<sampleID>.smooth.wig```<br>
Wig format files containing protein occupancy values at 10 base pairs across the whole genome.<br>
![pooled/<sampleID>.smooth.wig](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311139442.png)<br>


### Feature TSS coverage
```average_coverage.txt```<br>
"chr" column is the chromosome where each region is located.<br>
"start" column is the starting site of a region.<br>
"end" column is the ending site of a region.<br>
"coverage" column is average coverage value for each region.<br>
![average_coverage.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311175396.png)<br>

### Feature windowed protection score
```WPS.txt```<br>
"chr" column is the chromosome where each region is located.<br>
"start" column is the starting site of a region.<br>
"end" column is the ending site of a region.<br>
"long_WPS" column is windowed protection score of long reads for each region.<br>
"short_WPS" column is windowed protection score of short reads for each region.<br>
![WPS.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311508400.png)<br>

### Feature end motif frequency and diversity
```all_motifs_frequency.txt```<br>
Summary motif frequency for all the regions.<br>
![all_motifs_frequency.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311555448.png)<br>
<br>
```all_motifs_mds.txt```<br>
Summary motif diversity score for all the regions.<br>
![all_motifs_mds.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311576594.png)<br>
<br>
```region_<index>_motif_frequency_and_mds.txt```<br>
Motif frequency and diversity score for one region.<br>
![region_<index>_motif_frequency_and_mds.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311608092.png)<br>

### Feature fragmentation profile
```Fragmentation_Profile.txt```<br>
"seqnames" column is the chromosome where each region is located.<br>
"start" column is the starting site of a region.<br>
"end" column is the ending site of a region.<br>
"short" column is number of short reads identified in a region.<br>
"long" column is number of long reads identified in a region.<br>
"ratio" column is short/long ratio for each region.<br>
![Fragmentation_Profile.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311645495.png)<br>

### Feature nucleosome profile
```NucleosomeProfile.txt```<br>
"site_name" column is the name of each site list.<br>
"mean_coverage" column is the mean coverage for each site list.<br>
"central_coverage" column is the central coverage for each site list.<br>
"amplitude" column is the amplitude for each site list.<br>
![NucleosomeProfile.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311675968.png)<br>
<br>
```plots/<site_list>.pdf```<br>
Coverage profile for input bam file for given site lists.<br>


### Feature promoter fragmentation entropy
```PFE.txt```<br>
"TSS_ID" column is the id of each TSS 2k region.<br>
"PFE" column is promoter fragmentation entropy value for each TSS 2k region.<br>
![PFE.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311706499.png)<br>

### Feature orientation-aware cfDNA fragmentation
```OCF.txt```<br>
"TSS_ID" column is the id of each TSS 2k region.<br>
"OCF_Ratio" column is orientation-aware cfDNA fragmentation ratio for each TSS 2k region.<br>
![OCF.txt](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/1721311738958.png)<br>



## Versions of packages in our environment:
### Python:
```
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
```
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
```

### BiocManager packages
```
HMMcopy
GenomeInfoDb 
GenomicRanges
Rsamtools
GenomicAlignments
biovizBase
BSgenome.Hsapiens.UCSC.hg19
BSgenome.Hsapiens.UCSC.hg38
```
