# cfDNAanalyzer
cfDNAanalyzer (<ins>c</ins>ell-<ins>f</ins>ree <ins>DNA</ins> sequencing data <ins>analyzer</ins>) is a toolkit for cfDNA genomic sequencing data analysis which includes three core modules: <br>
(1) ***Feature Extraction***, which extracts multiple genomic and fragmentomic features at whole-genome or specific genomic-region levels;<br> 
(2) ***Feature Processing and Selection***, allowing for refinement and optimization of extracted features, suporting multiple processing and selection methods;<br> 
(3) ***Machine Learning Model Building***, supporting the development of both single-modality and multiple-modality predictive models for disease detection and classification.<br>
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
    <li><a href="#Structure-of-output-directory">Structure of output directory</a></li>
    <li><a href="#Features">Features</a></li>
    <li><a href="#Feature-Processing-and-Selection">Feature Processing and Selection</a></li>
    <li><a href="#Two-class-Machine-Learning">Two-class Machine Learning Models</a></li>
    <li><a href="#Multi-class-Machine-Learning">Multi-class Machine Learning Models</a></li>
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
Please ensure [<ins>samtools (v1.3.1)</ins>](https://github.com/samtools/samtools), [<ins>bedtools (v2.29.2)</ins>](https://bedtools.readthedocs.io/en/latest/index.html), and [<ins>deeptools (3.5.1)</ins>](https://github.com/deeptools/deepTools) are in your environment. Then, you can install the toolkit following the steps below ( R(>= 4.2.0) is required ):
```ruby
git clone https://github.com/LiymLab/cfDNAanalyzer.git
cd cfDNAanalyzer/
conda env create -f environment.yml
conda activate cfDNAanalyzer
Rscript install_R_packages.R
```
### Supported features
### 1.  Features for whole genome:
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1  <ins>C</ins>opy <ins>N</ins>umber <ins>A</ins>lterations (CNA) ([<ins>Adalsteinsson *et al, Nat. Commun.*, 2017</ins>](https://www.nature.com/articles/s41467-017-00965-y))<br>
* Copy number alterations comprise deletions or amplifications of a particular region of the genome, with a size as low as a few kilobases up to entire chromosomes. 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2  <ins>E</ins>nd <ins>M</ins>otif frequency and diversity (EM) ([<ins>Lee *et al, PNAS*, 2018</ins>](https://www.pnas.org/doi/abs/10.1073/pnas.1815031116) ; [<ins>Jiang *et al, Cancer Discov.*, 2020</ins>](https://doi.org/10.1158/2159-8290.CD-19-0622))
* End motifs are the terminal n-nucleotide sequence (4-mer end motif in this toolkit) at each 5′ fragment end of cfDNA molecules. End motif frequency refers to the frequency of all 256 4-mer end motifs.<br>
* End motif diversity is the normalized Shannon entropy of the categorical distribution of all possible 4-mer end-motifs of all cfDNA fragments.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.3  <ins>F</ins>ragmentation <ins>P</ins>rofile (FP) ([<ins>Cristiano *et al, Nature*, 2019</ins>](https://doi.org/10.1038/s41586-019-1272-6))
* Fragmentation profile describes fragmentation patterns of cfDNA across the genome, which is the fraction of short cfDNA fragments (100–150 bp) to long cfDNA fragments (151–220 bp) for each 5Mb window across the genome.

### 2.  Features for specific regions: 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.1  <ins>N</ins>ucleosome <ins>O</ins>ccupancy and <ins>F</ins>uzziness (NOF) ([<ins>Li *et al, Genome Med.*, 2024</ins>](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01280-6))
* Nucleosome occupancy reflects the frequency with which nucleosomes occupy a given genomic region in a cell population. Nucleosome occupancy for a specific region is calculated as the average occupancy values of all based in this region.<br>
* Nucleosome fuzziness is defined as the deviation of nucleosome positions within a region in a cell population and could reflect cell heterogeneity at the chromatin level. Nucleosome fuzziness for a specific region is calculated as the average fuzziness of all the nucleosomes whose center is located within the region.<br>
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.2 <ins>N</ins>ucleosome <ins>P</ins>rofile (NP) ([<ins>Doebley *et al, Nat. Commun.*, 2022</ins>](https://www.nature.com/articles/s41467-022-35076-w))
* Nucleosome profile is a composite coverage profile computed as the mean of the GC-corrected cfDNA fragment midpoint coverage across a set of sites (Binding sites sets of 377 transcription factors as the dafalut). For each set of sites, three features are identified from the composite coverage profile:
   * (1) The average coverage value from ± 30 bp of the central site of each set (central coverage).
   * (2) The average coverage value from ± 1000 bp of the central site of each set (mean coverage).
   * (3) The overall nucleosome peak amplitude is calculated by using a Fast Fourier Transform on the window ± 960 bp from the site and taking the amplitude of the 10th frequency term (amplitude).
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.3  <ins>W</ins>indowed <ins>P</ins>rotection <ins>S</ins>core (WPS) ([<ins>Snyder *et al, Cell*, 2016</ins>](https://doi.org/10.1016/j.cell.2015.11.050))
* Windowed protection score is the number of DNA fragments completely spanning a window centered at a given genomic coordinate minus the number of fragments with an endpoint within that same window. WPS can be calculated using long fragments (120–180 bp; 120 bp window) or short fragments (35–80 bp; 16 bp window). The WPS for a specific region is defined as the average WPS of all bases in this region. High WPS values indicate increased protection of DNA from digestion while low values indicate that DNA is unprotected.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.4  <ins>O</ins>rientation-aware <ins>C</ins>fDNA <ins>F</ins>ragmentation (OCF) ([<ins>Sun *et al, Genome Res.*, 2019</ins>](https://genome.cshlp.org/content/29/3/418.long))
* Orientation-aware cfDNA fragmentation is the differences of read densities of the upstream and downstream fragment ends in specific genomic regions.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.5  <ins>E</ins>nd <ins>M</ins>otif frequency and diversity for <ins>R</ins>egions (EMR) ([<ins>Serpas *et al, PNAS*, 2018</ins>](https://www.pnas.org/doi/abs/10.1073/pnas.1815031116) ; [<ins>Jiang *et al, Cancer Discov.*, 2020</ins>](https://doi.org/10.1158/2159-8290.CD-19-0622))
* We introduced end motif frequency and diversity for regions, which is defined as the frequency and diversity of all 256 4-mer end motifs for each region.<br>
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.6  <ins>F</ins>ragmentation <ins>P</ins>rofile for <ins>R</ins>egions (FPR) ([<ins>Cristiano *et al, Nature*, 2019</ins>](https://doi.org/10.1038/s41586-019-1272-6))
* We introduced fragmentation profile for regions, which is defined as the fraction of short cfDNA fragments (100–150 bp) to long cfDNA fragments (151–220 bp) for each region.

### 3.  Features for transcription start sites (TSSs): 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.1 <ins>P</ins>romoter <ins>F</ins>ragmentation <ins>E</ins>ntropy (PFE) ([<ins>Esfahani *et al, Nat. Biotechnol.*, 2022</ins>](https://doi.org/10.1038/s41587-022-01222-4))<br>
* Promoter fragmentation entropy quantifies the diversity of cfDNA fragment lengths around the TSSs of genes. It is calculated by a modified Shannon index for lengths of cfDNA fragment where both ends fell within ±1 kb of the TSS. Then this cfDNA entropy measure is adjusted using a Dirichlet-multinomial model for normalization.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2 <ins>TSS</ins> <ins>C</ins>overage (TSSC) ([<ins>Ulz *et al, Nat. Genet.*, 2016</ins>](https://www.nature.com/articles/ng.3648))
* TSS coverage refers to the cfDNA sequencing read coverage in the regions surrounding TSSs.

### Supported feature processing methods and machine learning models

### 1. Feature processing 

#### 1.1 Missing data removal
* In this step, columns containing missing values are removed to ensure the data is ready for downstream analysis.

#### 1.2 Data standardization

* This process transforms the data using one the following data standardization/normalization method, which is important for many machine-learning techniques. 

| Methods | Alias | Brief Introduction |
| :---: | :---: | :---: | 
| Z-Score | *Zscore* | Z-Score scales the features to mean = 0, standard deviation = 1. |
| Min-Max Scaling | *MMS* | Min-Max Scaling transforms features to a fixed range (typically [0, 1]) by subtracting the minimum value and dividing by the range (max - min). |  
| Robust Scaling | *RS* | Robust Scaling standardizes data by subtracting the median and dividing by the interquartile range (IQR), making it resistant to outliers. |
| Quantile Transformation | *QT* | Quantile Transformation non-linearly maps data to a target distribution (e.g., Gaussian or uniform) using quantiles, making the transformed data robust to outliers and strictly follow the desired distribution. |

### 2. Feature selection

Feature selection is an important step in machine leanring. It involves removing irrelevant or redundant features based on importance rankings, which can help reduce model complexity and enhance the model's ability to make accurate predictions. Here we provide **four methods (filter, embedded, wrapper, and hybrid methods)** for feature selection ( [<ins>Pudjihartono *et al, Front. Bioinform.*, 2022</ins>](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9580915/)), as follows:

#### 2.1 Filter methods

Filter methods apply statistical tests to assess and rank features for their relevance to the target variable, independent of the classification algorithm. These methods can be executed by computing the correlation between features and the target variable, subsequently eliminating features with low scores.

| Methods | Alias | Brief Introduction | Source |
| :---: | :---: | :---: | :---:                                                   |
| Information Gain | *IG* | Information Gain measures the reduction in uncertainty of the target variable when a specific feature is known. Higher information gain indicates a more informative feature for predicting the target variable. | [scikit-learn](https://scikit-learn.org/stable/index.html)                                                 |
| Chi-square Test | *CHI* | The Chi-square Test evaluates the independence between categorical features and the target variable. Features that show significant association with the target variable are selected based on their Chi-square statistics. | [scikit-learn](https://scikit-learn.org/stable/index.html)                                                 |
| Fisher's Score | *FS* | Fisher’s Score ranks features according to their ability to distinguish between classes by comparing the mean differences across classes relative to their variances. | [scikit-feature](https://github.com/jundongl/scikit-feature/tree/master) |
| FCBF Fast Correlation Filter | *FCBF* | The FCBF method quickly identifies relevant features by analyzing their correlation with the target variable while removing redundant features that contribute little additional information. | [FCBF_module](https://github.com/SantiagoEG/FCBF_module)     |
| Permutation Importance | *PI* | Permutation Importance assesses the significance of each feature by measuring the drop in model performance when the feature’s values are randomly permuted. Features causing significant degradation in performance are deemed important. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master)                                                  |
| Correlation Coefficient | *CC* | Correlation Coefficient quantifies the linear relationship between features and the target variable. Features with higher correlation coefficients are considered more predictive and are selected for modeling. | [pandas](https://pandas.pydata.org/)                                                       |
| Low Variance Filter | *LVF* | The Low Variance Filter eliminates features with minimal variability across samples, as such features are unlikely to be useful for distinguishing between classes. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master)                                                 |
| Mean Absolute Difference | *MAD* | Mean Absolute Difference evaluates the average difference in feature values between classes. Features with larger differences are prioritized for greater discriminatory power. | [pandas](https://pandas.pydata.org/)                                                       |
| Dispersion Ratio | *DR* | Dispersion Ratio assesses the spread or variability of a feature across samples. Features with distinct dispersion patterns across classes are selected as they are likely to contribute to better classification. | [numpy](https://numpy.org/)                                                        |
| Mutual Information | *MI* | Mutual Information measures the amount of information a feature shares with the target variable. Features with higher mutual information values are more relevant for model inclusion. | [scikit-learn](https://scikit-learn.org/stable/index.html)                                                 |
| ReliefF | *RLF* | ReliefF identifies features that consistently differentiate between instances of different classes by considering the nearest neighbors of each instance. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
| SURF | *SURF* | The SURF variant extends ReliefF by considering all instances within a predefined distance, enhancing its ability to capture more relevant features. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
| MultiSURF | *MSURF* | MultiSURF dynamically adjusts the neighborhood size during selection, allowing it to capture complex, multi-feature interactions. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
| TuRF | *TRF* | TuRF iteratively removes the least informative features, refining the feature set to focus on those with the highest discriminative power. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
#### 2.2 Embedded methods

Embedded methods incorporate feature selection directly into model training, allowing the learning algorithm to determine the most crucial features. Models such as Random Forest and LASSO regression effectively execute these methods, assigning importance to features throughout the learning phase.

| Methods | Alias | Brief Introduction | Source |
| :---: | :---: | :---: | :---: |
| Lasso | *LASSO* | Lasso performs feature selection by applying L1 regularization, which encourages sparsity in the model by driving the coefficients of less important features to zero. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| Ridge | *RIDGE* | Ridge applies L2 regularization, which penalizes large coefficients, thereby reducing the impact of less important features while retaining all features in the model. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| ElasticNet | *ELASTICNET* | ElasticNet combines the strengths of both Lasso (L1 regularization) and Ridge (L2 regularization), allowing for balanced feature selection with both sparsity and regularization. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| Random Forest Importance | *RF* | Random Forest Importance ranks features based on their contribution to model accuracy. The method leverages the structure of decision trees to identify the most predictive features during the training process. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |

#### 2.3 Wrapper methods

Wrapper methods choose features according to their impact on a selected classifier’s performance, systematically seeking the optimal subset of features. This is done by evaluating various feature combinations and identifying the subset that enhances the model’s performance, typically using techniques such as recursive feature elimination.

| Methods | Alias | Brief Introduction | Source |
| :---: | :---: | :---: | :---: |
| Forward Feature Selection | *FFS* | Forward Feature Selection starts with an empty model, sequentially adding features that improve model performance until no further gain is achieved. This method builds the optimal feature set step by step. | [mlxtend](https://rasbt.github.io/mlxtend/) |
| Backward Feature Selection | *BFS* | Backward Feature Selection begins with all features included and iteratively removes the least significant ones. The process continues until the removal of further features results in performance degradation. | [mlxtend](https://rasbt.github.io/mlxtend/) |
| Exhaustive Feature Selection | *EFS* | Exhaustive Feature Selection evaluates all possible feature combinations, ensuring the identification of the best subset. Though computationally intensive, it provides the most thorough feature selection. | [mlxtend](https://rasbt.github.io/mlxtend/)  |
| Recursive Feature Selection | *RFS* | Recursive Feature Selection recursively eliminates the least important features, refining the model by retaining only those that contribute the most to prediction accuracy. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| Boruta Algorithm | *BOR* | The Boruta Algorithm iteratively compares the importance of real features with shuffled versions, selecting all features that consistently prove to be more informative than random noise. | [boruta](https://www.jstatsoft.org/article/view/v036i11) |

#### 2.4 Hybrid methods

Hybrid approaches combine filter with wrapper or embedded methods to leverage their unique advantages. For instance, features are initially filtered using statistical tests, followed by a wrapper method that enhances the selection by examining interactions with the classifier.

| Methods | Alias | Brief Introduction |
| :---: | :---: | :---: |
| Filter-Wrapper | *FW* | Filter methods are firstly used to reduce the feature set based on relevance, and then wrapper methods refine the selection to identify the optimal features for the model. |
| Filter-Embedded | *FE* | Filter methods are firstly used to reduce the feature set based on relevance, and then embedded methods refine the selection to identify the optimal features for the model. |

### 3. Single modality machine learning model

In this section, we will use supervised machine learning to analyze the **individual feature type** of all samples. The classifiers in the toolkit include **Random Forest, K-Nearest Neighbors, Gaussian Naive Bayes, Logistic Regression, Support Vector Machine, and eXtreme Gradient Boosting**. The goal of this approach is to assess how effective each classifier is in accurately identifying disease status based on individual feature type. Our methods are applicable to both **binary and multi-class classification** scenarios.

### 4. Multiple modality machine learning model

This section aims to use **multiple feature types** to improve the model performance. Our toolkit include **Concatenation-based, Model-based, and Transformation-based methods** ([Reel *et al, Biotechnol. Adv.*, 2021](https://www.sciencedirect.com/science/article/abs/pii/S0734975021000458?via%3Dihub/)), applicable to both **binary and multi-class classification** tasks, to maximize model effectiveness across different scenarios. 

#### 4.1 Concatenation-based methods

Concatenation-based methods combine features from different sources by directly appending them into a single feature vector. This straightforward approach allows models to access all available data but relies heavily on the subsequent learning algorithm to handle the potentially high-dimensional space effectively.

#### 4.2 Model-based methods

Model-based integration methods create multiple intermediate models for the different omics data and then build a final model from various intermediate models ([Reel *et al, Biotechnol. Adv.*, 2021](https://www.sciencedirect.com/science/article/abs/pii/S0734975021000458?via%3Dihub/)). It boosts predictive accuracy by utilizing the strengths of individual models, enabling a deeper understanding of the interactions across different layers.

| Methods | Brief Introduction | 
| :---: | :---: | 
| Average Voting | Average Voting combines predictions from multiple models by calculating the mean of their outputs. This method can balance diverse feature impacts, reducing the effect of anomalies in individual features. | 
| Weighted Voting | Weighted Voting assigns weights to individual model predictions based on their performance, using either their proportion of the total AUC or inversely to their error rates. This strategy prioritizes more reliable features when integrating multiple feature sets, enhancing overall prediction accuracy. | 
| Majority Voting | Majority Voting selects the most common outcome among various model predictions. It ensures robustness by reducing the influence of sporadically erroneous model outputs. | 
| Stack Ensemble | Stack Ensemble uses a meta-model to learn how to best combine the predictions of multiple base models. This method effectively integrates diverse feature-derived predictions to optimize final decision-making processes. | 

#### 4.3 Transformation-based methods

Transformation-based methods transform each of the omics datasets firstly into graphs or kernel matrices and then combine all of them into one before constructing a model ([Reel *et al, Biotechnol. Adv.*, 2021](https://www.sciencedirect.com/science/article/abs/pii/S0734975021000458?via%3Dihub/)). This approach supports the integration of varied data by standardizing how they are represented.

| Category | Methods | Brief Introduction | Reference |
| :---: | :---: | :---: | :---: |
| Dimension reduction | PCA | PCA reduces dimensionality by transforming features into a set of linearly uncorrelated components. It simplifies multi-feature integration by focusing on components that explain the most variance, thereby enhancing model efficiency and clarity. | [Subramanian *et al, Bioinform. Biol. Insights*, 2020](https://journals.sagepub.com/doi/10.1177/1177932219899051/) |
| Kernel matrix | Linear Kernel | The Linear Kernel measures direct linear relationships between features. It aids in multi-feature integration by emphasizing linear associations, making it suitable for linearly separable data. |  |
| Kernel matrix | Polynomial Kernel | The Polynomial Kernel allows the capture of interactions between features at different degrees of complexity. It enhances multi-feature integration by providing a flexible framework to model nonlinear relationships in the data. |  |
| Kernel matrix | Radial Basis Function (RBF) Kernel | The RBF Kernel maps features into a higher-dimensional space using a radial basis function, enabling effective classification of non-linearly separable datasets. This kernel is crucial for multi-feature integration as it can handle complex and non-linear interactions between features. |  |
| Kernel matrix | Sigmoid Kernel | The Sigmoid Kernel projects data using a sigmoid function, similar to neural network activation functions. It supports multi-feature integration by transforming features into formats that highlight threshold-based classifications. |  |
| Network | Similarity network fusion (SNF) | SNF integrates multiple types of data by fusing similarity networks, reinforcing common structural features. It excels in multi-feature integration by constructing a holistic network view, revealing deep insights across combined datasets. | [Wang *et al, Nat. Methods*, 2014](https://www.nature.com/articles/nmeth.2810) |

### Usage
```ruby
bash cfDNAanalyzer.sh -I <InputFile> -o <OutputDirectory> -F <Features> [Options]
```
#### Options: 
```
-------------------- Options for feature extraction --------------------

-----General options -----
  -I  FILE      A text file containing all input BAM files with one BAM file per line. 
                BAM files generated using both Bowtie2 and BWA are accepted.  
  -o  DIR       Output directory for all the results. Default: [./]
  -F  STR       Features to extract, including CNA, NOF, WPS, EM, EMR, FP, FPR, NP, OCF, PFE, and TSSC.
                Features should be set as a string separated by comma, e.g., CNA,NOF. 
                Default: All available features will be extracted.
                The detailed description of each feature can be accessed at https://github.com/LiymLab/cfDNAanalyzer. 
                Note: The following features are specifically designed for paired-end sequencing data: FP, FPR, EM, EMR, NP, PFE, and OCF.
  -g  STR       Genome version of input BAM files (hg19/hg38). Default: [hg38] 
  -b  FILE      A BED3 file specifying the regions to extract features.
                The file should contain three TAB-delimited columns: chromosome start end.
  -s  STR       Sequencing method of input BAM files (single/pair). Default: [pair]
  -t  INT       Number of threads to use. Default: [1]              

----- Options specific for Copy Number Alterations (CNA) -----
  -B     INT    Bin size in kilobases (10, 50, 500, or 1000). Default: [1000]
  --CNA  STR    Additional parameter setting for CNA analysis. 
                The full parameter list is available by running: Rscript cfDNAanalyzer/ichorCNA/ichorCNA/scripts/runIchorCNA.R --help. [optional]

----- Options specific for Nucleosome Occupancy and Fuzziness (NOF) -----
  --NOF  STR    Additional parameter setting for NOF analysis. 
                The full parameter list is available by running: python cfDNAanalyzer/DANPOS3/danpos.py dpos -h. [optional]

----- Options specific for Windowed Protection Score (WPS) -----
  -x  INT       Min fragment length used for long fragments WPS calculation. Default: [120]
  -X  INT       Max fragment length used for long fragments WPS calculation. Default: [180]
  -w  INT       Window size for long fragments WPS calculation. Default: [120]
  -m  INT       Min fragment length used for short fragments WPS calculation. Default: [35]
  -M  INT       Max fragment length used for short fragments WPS calculation. Default: [80]
  -W  INT       Window size for short fragments WPS calculation. Default: [16]

----- Options specific for End Motif frequency and diversity (EM) -----
  -f  FILE      Reference genome in FASTA format. 
                For example, hg38 reference genome can be downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz.

----- Options specific for Nucleosome Profile (NP) -----
  -l  DIR       Directory containing a list of files with each file for a set of sites.
                The file must have at least two columns with the following column names: 
                "Chrom" for the chromosome name and "position" for the site positions. 
                If not provided, the 377 TF binding site lists (cfDNAanalyzer/Griffin/Ref_hg19/sites or cfDNAanalyzer/Griffin/Ref_hg38/sites) from the original Nucleosome Profile paper will be used . 

----- Options for Promoter Fragmentation Entropy (PFE) -----
  -T     FILE   A TAB-delimited TSS information file without any header.
                The file must have six columns: (1) chromosome, (2) 1-base TSS coordinate, (3) gene name, (4) strand (+1/-1) and (5) TSS ID (e.g., for genes with multiple TSS, this could be geneName_1, geneName_2, etc.)
                If not provided, TSSs (cfDNAanalyzer/Epic-seq/code/priordata/sample_hg19.txt or cfDNAanalyzer/Epic-seq/code/priordata/sample_hg38.txt) from the original Promoter Fragmentation Entropy paper will be used.
  --PFE  STR    Addtional parameter setting for PFE analysis. 
                The full parameter list is available by running: Rscript cfDNAanalyzer/Epic-seq/code/runEPIC.R -h.[optional]

----- Options for TSS Coverage (TSSC) -----
  -u                    INT   Number of base pairs upstream of TSS used for calculating TSSC. Default: [1000]
  -d                    INT   Number of base pairs downstream of TSS used for calculating TSSC. Default: [1000]
  -S                    FILE  A BED6 file specifying the coordinates of TSSs used for calculating TSSC. 
                              If not provided, all TSSs of refSeq genes will be used.
  -n                    STR   Methods to normalize the number of reads per bin for the "bamCoverage" command of deeptools (RPKM,CPM,BPM,RPGC) . Default: [RPKM]
  --bamCoverage         STR   Additional parameter seting for the "bamCoverage" command of deeptools.
                              The full parameter list is available by running: bamCoverage -h. [optional] 
  --multiBigwigSummary  STR   Additional parameter seting for the "multiBigwigSummary" command of deeptools.
                              The full parameter list is available by running: multiBigwigSummary -h.[optional].
                          

-------------------- Options for downstream analysis --------------------
  --noDA                LOGIC  Skip all the downstream analysis.
  --noML                LOGIC  Skip machine learning model building, only feature processing and selection will be conducted.
                               Methods in feature selection will leverage the label information of all samples.
  --labelFile           FILE   Label information file for all the samples.
                               This file must have two columns. The sample column contains the basename of all the samples. The label column contains all the samples' label information (like 0,1 for two class and 0,1,2 for three class).

----- Options for feature processing and selection -----
  --standardM                    STR    The standardization method (Zscore, MMS, RS, or QT). Default:   [Zscore]
                                        The detailed description of each method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --filterMethod                 STR    Filter methods employed for feature selection (IG CHI FS FCBF PI CC LVF MAD DR MI RLF SURF MSURF TRF).
                                        Methods should be set as a string separated by space, e.g., MR IG CHI. 
                                        Default: [IG].
                                        The detailed description of each filter method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --filterFrac                   FLOAT  Fraction of features to retain when employing the filter method. Default: [0.2]
  --wrapperMethod                STR    Wrapper methods employed for feature selection (FFS BFS EFS RFS BOR).
                                        Methods should be set as a string separated by space, e.g., BOR RFS. 
                                        Default: [BOR].
                                        The detailed description of each wrapper method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --wrapperFrac                  FLOAT  Fraction of features to retain when employing the wrapper method. Default: [0.2]
  --embeddedMethod               STR    Embedded methods employed for feature selection (LASSO RIDGE ELASTICNET RF).
                                        Methods should be set as a string separated by space, e.g., LASSO RIDGE. 
                                        Default: [LASSO]
                                        The detailed description of each embedded method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --embeddedFrac                 FLOAT  Fraction of features to retain when employing the embedded method. Default: [0.2]
  --hybridType                   STR    Two methods used for hybrid method (FW/FE). Default: [FE]
  --hybridMethod1,hybridMethod2  STR    Subtype methods designated for method 1 and method 2 in the "--hybridType".
                                        Methods should be set as a string separated by space, e.g., BOR RFS.
                                        Default: [IG,BOR] for FW; [IG,LASSO] for FE.
  --hybridFrac1,hybridFrac2       FLOAT Fraction of features to retain for method 1 and method 2 in the "--hybridType". Default: [0.2,0.2]
  
----- Options for machine learning model building -----
  --classNum          INT   Number of classification categories (2 or more). Default: [2]
  --cvSingle          STR   Cross-validation method applied in single modality, options include leave-one-out (LOO) or K-fold (KFold). Default: [LOO]
  --nsplitSingle      INT   Number of folds designated for K-fold cross-validation in single modality. Default: [5]
  --classifierSingle  STR   Classifiers employed in single modality (KNN SVM RandomForest GaussianNB LogisticRegression XGB).
                            Classifiers should be set as a string separated by space, e.g., KNN XGB.
                            Default: All available classifiers will be applied.
  --cvMulti           STR   Cross-validation method applied in multiple modalities, options include leave-one-out (LOO) or K-fold (KFold). Default: [LOO]
  --nsplitMulti       INT   Number of folds designated for K-fold cross-validation in multiple modalities. Default: [5]
  --classifierMulti   STR   Classifiers employed in multiple modalities (KNN SVM RandomForest GaussianNB LogisticRegression XGB).
                            Classifiers should be set as a string separated by space, e.g., KNN XGB.
                            Default: All available classifiers will be applied.  
```

### Run the usage example
You can directly run cfDNAanalyzer by ```cfDNAanalyzer.sh``` provided in the ```cfDNAanalyzer/``` directory with the example files available at Zenodo (DOI: ([<ins>10.5281/zenodo.13369741</ins>](https://doi.org/10.5281/zenodo.13369741))).  

```ruby
bash cfDNAanalyzer.sh -I ./example/bam_input.txt -F CNA,FP,OCF -g hg19 -b ./example/test.bed -f ./example/example.fa --filterMethod MR IG CHI --wrapperMethod BOR --embeddedMethod LASSO RIDGE --classNum 2 --cvSingle LOO --classifierSingle KNN SVM RandomForest --cvMulti LOO --classifierMulti KNN SVM RandomForest > ./cfDNAanalyzer.log
```

## Output files
### Structure of output directory
```shell
├── Features
│   ├── CNA.csv
│   ├── EM_motifs_frequency.csv
│   ├── EM_motifs_mds.csv
│   ├── FP_fragmentation_profile.csv
│   ├── NOF_meanfuziness.csv
│   ├── NOF_occupancy.csv
│   ├── NP_amplitude.csv
│   ├── NP_central_coverage.csv
│   ├── NP_mean_coverage.csv
│   ├── WPS_long.csv
│   ├── WPS_short.csv
│   ├── OCF.csv
│   ├── EMR_aggregated_motif_frequency.csv
│   ├── EMR_aggregated_mds.csv
│   ├── EMR_region_motif_frequency.csv
│   ├── EMR_region_mds.csv
│   ├── FPR_fragmentation_profile_regions.csv
│   ├── PFE.csv
│   └── TSSC_average_coverage.csv
├── Feature_Processing_and_Selection
│   ├── Feature_Processing
│   │   └── [FeatureName].csv
│   └── Feature_Selection
│       ├── [FeatureName]_[embeddedMethod]_selectd.csv
│       ├── [FeatureName]_[filterMethod]_selectd.csv
│       └── [FeatureName]_[wrapperMethod]_selectd.csv
└── Machine_Learning
    ├── single_modality
    │   ├── [FeatureName1]
    │   │   ├── single_modality_metrics.csv
    │   │   └── single_modality_probabilities.csv
    │   │   ......
    │   ├── [FeatureName2]
    │   └── [FeatureNameN]
    └── multiple_modality
        ├── Concatenation_based
        │   ├── multiple_modality_metrics.csv
        │   └── multiple_modality_probabilities.csv
        ├── Model_based
        │   ├── multiple_modality_metrics.csv
        │   └── multiple_modality_probabilities.csv
        └── Transformation_based
            ├── multiple_modality_metrics.csv
            └── multiple_modality_probabilities.csv
```
### Features
Output files under `Features` directory consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification.
#### Copy Number Alterations (CNA)
In ```CNA.csv```, columns after `label` contain the estimated copy number for each bin. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each bin.
```r
sample,label,chr10_100000001_100100000,chr10_10000001_10100000,chr10_1000001_1100000,...
sample1,1,2,2,2,...
sample2,1,2,2,2,...
sample3,0,2,2,2,...
sample4,0,3,3,3,...
```

#### End Motif frequency and diversity of whole genome (EM)
In ```EM_motifs_frequency.csv```, columns after `label` contain motif frequency for input BAM file across the whole genome. Column name `[Motifs]` of these columns specifies different end motifs.
```r
sample,label,AAAA,AAAC,AAAG,...
sample1,1,0.010256,0.003589,0.005923,...
sample2,1,0.009953,0.003414,0.005793,...
sample3,0,0.010264,0.003594,0.005878,...
sample4,0,0.00947,0.00344,0.005482,...
```
In ```EM_motifs_mds.csv```, column `MDS` after `label` contain motif frequency diversity scores for input BAM file across the whole genome.
```r
sample,label,MDS
sample1,1,0.965699
sample2,1,0.967610
sample3,0,0.969939
sample4,0,0.974932
```

#### Fragmentation Profile of whole genome (FP)
In ```Fragmentation_Profile.txt```, columns after `label` contain ratio of short to long fragments across the whole genome with a window of 100kb. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each window.
```r
sample,label,chr10_100000001_100100000,chr10_10000001_10100000,chr10_1000001_1100000,...
sample1,1,0.143288,0.182047,0.205779,...
sample2,1,0.232316,0.243959,0.238103,...
sample3,0,0.167793,0.16843,0.187116,...
sample4,0,0.131223,0.180202,0.200336,...
```

#### Nucleosome Occupancy and Fuzziness (NOF)
In ```NOF_meanfuziness.csv``` and ```NOF_occupancy.csv```, columns after `label` contain average fuzziness value (```NOF_meanfuziness.csv```) or occupany value (```NOF_occupancy.csv```) for each region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,5.10545,15.2444,13.1614,...
sample2,1,14.3248,15.8161,21.7031,...
sample3,0,8.85557,15.4003,12.8956,...
sample4,0,91.5162,132.746,135.676,...
```
#### Nucleosome Profile (NP)
In ```NP_amplitude.csv```, ```NP_central_coverage.csv``` and ```NP_mean_coverage.csv```, columns after `label` contain nucleosome peak amplitude (```NP_amplitude.csv```), central coverage (```NP_central_coverage.csv```), and the mean coverage (```NP_mean_coverage.csv```) for each site set in the input directory. Column name `[site].10000.txt` of these columns specifies different sites set.
```r
sample,label,AHR.10000.txt,AR.10000.txt,ARID3A.10000.txt,...
sample1,1,0.98707,1.00713,0.97745,...
sample2,1,0.99656,1.01126,0.98364,...
sample3,0,0.97693,1.00172,0.97278,...
sample4,0,0.98448,1.01047,0.98699,...
```
#### Windowed Protection Score (WPS)
In ```WPS_long.csv``` and ```WPS_short.csv```, columns after `label` contain WPS for long fragments (```WPS_long.csv```), WPS for short fragments (```WPS_short.csv```) for each region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,-0.453047,-1.16434,-0.827173,...
sample2,1,-1.12138,-0.908591,-1.54246,...
sample3,0,-0.934066,-0.972028,-0.815684,...
sample4,0,-7.13536,-8.07143,-7.43057,...
```

#### Orientation-aware CfDNA Fragmentation (OCF)
In ```OCF.csv```, columns after `label` contain OCF values for each region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,0,434.782608,0,...
sample2,1,0,0,44.642857,...
sample3,0,0,0,0,...
sample4,0,0,-270.564042,52.910052,...
```

#### End Motif frequency and diversity for regions (EMR)
Our toolkit outputs two types of EMR features at different levels:<br>
(1) Motif frequency and diversity for all input regions aggregated together.<br>
```EMR_aggregated_motif_frequency.csv```<br>
```r
sample,label,AAAA,AAAC,AAAG,...
sample1,1,0.010256,0.003589,0.005923,...
sample2,1,0.009953,0.003414,0.005793,...
sample3,0,0.010264,0.003594,0.005878,...
sample4,0,0.00947,0.00344,0.005482,...
```
```EMR_aggregated_mds.csv```<br>
```r
sample,label,MDS
sample1,1,0.965699
sample2,1,0.967610
sample3,0,0.969939
sample4,0,0.974932
```
&nbsp;<br>
(2) Motif frequency and diversity for each region separately.<br>
```EMR_region_motif_frequency.csv```<br>
```r
sample,label,chr10_100026951_100028952_AAAA,chr10_100026951_100028952_AAAC,chr10_100026951_100028952_AAAG,...
sample1,1,0.000000,0.000983,0.000983,...
sample2,1,0.001222,0.000815,0.000407,...
sample3,0,0.001343,0.003358,0.002686,...
sample4,0,0.000595,0.001123,0.002114,...
```
```EMR_region_mds.csv```<br>
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,0.905673,0.951263,0.948304,...
sample2,1,0.928867,0.941905,0.948625,...
sample3,0,0.932926,0.946175,0.950612,...
sample4,0,0.937257,0.951996,0.954475,...
```

#### Fragmentation Profile for regions (FPR)
In ```FPR_fragmentation_profile_regions.csv```, columns after `label` contain ratio of short to long fragments for each input region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100190117_100192118,...
sample1,1,0.622159,0.205120,0.255837,...
sample2,1,0.856960,1.246473,0.596067,...
sample3,0,2.475976,0.387795,0.469562,...
sample4,0,1.254673,0.769168,0.870592,...
```


#### Promoter Fragmentation Entropy (PFE)
In ```PFE.csv```, columns after `label` contain PFE value for the 2kb surrounding TSS. Column name `[TSS]` of these columns specifies different TSS ID.
```r
sample,label,A1BG_1,A1CF_1,A2ML1_1,...
sample1,1,0.342051,0.342401,0.345179,...
sample2,1,0.378536,0.378173,0.376558,...
sample3,0,0.382465,0.380891,0.377782,...
sample4,0,0.365,0.365356,0.366375,...
```

#### TSS Coverage (TSSC)
In ```TSSC_average_coverage.csv```, columns after `label` contain cfDNA fragment coverage value for each TSS surrounding regions. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each TSS surrounding regions.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,1.642751,1.212495,1.000397,...
sample2,1,1.406582,1.938656,1.691248,...
sample3,0,1.220205,1.239343,2.053586,...
sample4,0,1.180354,1.51295,0.361517,...
```

### Feature Processing and Selection
Output files under `Feature_Processing_and_Selection/Feature_Processing` and `Feature_Processing_and_Selection/Feature_Selection` directory consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification.

#### Feature Processing

In `[FeatureName].csv`, columns after `label` contain **processed** feature data for each sample.
```r
sample,label,feature1,feature2,feature3,...
sample1,1,-0.710433,-0.206147,-0.419211,...
sample2,1,0.354042,1.435867,0.320478,...
sample3,0,-0.417439,-0.567277,-0.846278,...
sample4,0,-0.85468,-0.255078,-0.54376,...
```

#### Feature Selection

In`[FeatureName]_[embeddedMethod]_selectd.csv`,`[FeatureName]_[filterMethod]_selectd.csv` and `[FeatureName]_[wrapperMethod]_selectd.csv`, columns after `label` contain **selected** feature data for each sample.

```r
sample,label,feature1,feature2,feature3,...
sample1,1,0.612106,0.9234,0.912867,...
sample2,1,0.002909,-0.722122,-0.024669,...
sample3,0,0.626594,0.970688,0.588356,...
sample4,0,-0.967271,-0.476624,-2.259162,...
```

### Two-class Machine Learning

#### Single Modality/[FeatureName]

`single_modality_metrics.csv` contains the classifier and feature selection method used, followed by various performance metrics such as accuracy, precision, recall, F1 score, AUC (Area Under the Curve), computation time and memory usage (peak memory).

```r
FS_Combination,Classifier,accuracy,precision,recall,f1,auc,TotalTime_sec,PeakMemory_MB
filter_IG_0.023618328,KNN,0.48,0.481481,0.52,0.5,0.4368,919.7527,409.4531
filter_IG_0.023618328,SVM,0.5,0.5,0.56,0.528302,0.4768,1835.6222,401.8867
filter_CHI_0.023618328,KNN,0.6,0.6,0.6,0.6,0.6464,42.5313,485.0039
filter_CHI_0.023618328,SVM,0.64,0.612903,0.76,0.678571,0.4528,33.3784,559.9922
filter_FS_0.023618328,KNN,0.5,0.0,0.0,0.0,0.5,3.4107,601.1992
```

`single_modality_probabilities.csv` contains samples IDs and their true labels, followed by the classifier, feature selection method used and predicted probabilities for each class.

```r
SampleID,TrueLabel,FS_Combination,Classifier,Prob_Class0,Prob_Class1
sample1,0,KNN,filter_IG_0.023618328,0.6,0.4
sample2,1,KNN,filter_IG_0.023618328,0.4,0.6
sample3,1,KNN,filter_IG_0.023618328,0.8,0.2
sample4,0,KNN,filter_IG_0.023618328,0.4,0.6
sample5,1,KNN,filter_IG_0.023618328,0.6,0.4
```

#### Multiple Modality/Concatenation based, Multiple Modality/Model based, Multiple Modality/Transformation based
`multiple_modality_metrics.csv` contains the fusion type, fusion method, classifier and feature selection method used, followed by various performance metrics such as accuracy, precision, recall, F1 score, AUC (Area Under the Curve), computation time and memory usage (peak memory).

```r
FusionType,FusionMethod,FS_Combination,Classifier,accuracy,precision,recall,f1,auc,TotalTime_sec,PeakMemory_MB
concat,concat,filter_DR_0.2,KNN,0.52,0.52381,0.44,0.478261,0.5616,0.0002,232.9531
concat,concat,filter_DR_0.2,SVM,0.4,0.4,0.4,0.4,0.3776,0.0001,237.3281
concat,concat,filter_DR_0.2,KNN,0.52,0.52381,0.44,0.478261,0.5616,0.0002,232.9531
concat,concat,filter_DR_0.2,SVM,0.4,0.4,0.4,0.4,0.3776,0.0001,237.3281
```

`multiple_modality_probabilities.csv` contains samples IDs and their true labels, followed by the  fusion type, fusion method, classifier, feature selection method used and predicted probabilities for each class.

```r
SampleID,TrueLabel,FusionType,FusionMethod,FS_Combination,Classifier,Prob_Class0,Prob_Class1
sample1,1,concat,concat,filter_DR_0.2,KNN,0.4,0.6
sample2,1,concat,concat,filter_DR_0.2,KNN,0.6,0.4
sample3,1,concat,concat,filter_DR_0.2,KNN,0.4,0.6
sample4,1,concat,concat,filter_DR_0.2,KNN,0.6,0.4
sample5,1,concat,filter_DR_0.2,KNN,concat,0.4,0.6
```

### Multi-class Machine Learning

#### Single modality
`single_modality_metrics.csv` contains classifier and feature selection method used, followed by various performance metrics such as accuracy, macro-precision, macro-recall, macro-f1 score, computation time and memory usage (peak memory).
```r
FS_Combination,Classifier,accuracy,precision_macro,recall_macro,f1_macro,TotalTime_sec,PeakMemory_MB
wrapper_BOR_0.021070375,KNN,0.473684,0.157895,0.333333,0.214286,509.695,992.7558
wrapper_BOR_0.021070375,SVM,0.421053,0.488889,0.342593,0.327778,421.3054,991.1528
wrapper_BOR_0.021070375,RandomForest,0.578947,0.427778,0.546296,0.472222,432.2937,991.2022
wrapper_BOR_0.021070375,GaussianNB,0.421053,0.305556,0.527778,0.385185
wrapper_BOR_0.021070375,LogisticRegression,580.3495,990.9873,0.473684,0.488889,0.425926,0.416667,551.8989,992.4312
wrapper_BOR_0.021070375,XGB,0.421053,0.3125,0.342593,0.297778,571.5343,963.2724
```

`single_modality_probabilities.csv` contains samples IDs and their true labels, followed by the classifier, feature selection method used and predicted probabilities for each class.

```r
SampleID,TrueLabel,FS_Combination,Classifier,Prob_Class0,Prob_Class1,Prob_Class2
sample1,wrapper_BOR_0.021070375,0,KNN,0.2,0.8,0.0
sample2,wrapper_BOR_0.021070375,2,KNN,0.2,0.6,0.2
sample3,wrapper_BOR_0.021070375,2,KNN,0.0,0.6,0.4
sample4,wrapper_BOR_0.021070375,0,KNN,0.2,0.8,0.0
```

#### Multiple Modality/Concatenation based, Multiple Modality/Model based, Multiple Modality/Transformation based
`multiple_modality_metrics.csv` contains the fusion type, fusion method, classifier and feature selection method used, followed by various performance metrics such as accuracy, macro-precision, macro-recall, macro-f1 score, computation time and memory usage (peak memory).

```r
FusionType,FusionMethod,FS_Combination,Classifier,accuracy,precision_macro,recall_macro,f1_macro,TotalTime_sec,PeakMemory_MB
concat,concat,filter_DR_0.2,KNN,0.52,0.52381,0.44,0.478261,0.5616,232.9531
concat,concat,filter_DR_0.2,SVM,0.4,0.4,0.4,0.4,0.0001,237.3281
concat,concat,filter_DR_0.2,KNN,0.52,0.52381,0.44,0.5616,0.0002,232.9531
concat,concat,filter_DR_0.2,SVM,0.4,0.4,0.4,0.4,0.0001,237.3281
```

`multiple_modality_probabilities.csv` contains samples IDs and their true labels, followed by the  fusion type, fusion method, classifier, feature selection method used and predicted probabilities for each class.

```r
SampleID,TrueLabel,FusionType,FusionMethod,FS_Combination,Classifier,Prob_Class0,Prob_Class1,Prob_Class2
sample1,1,concat,concat,filter_DR_0.2,KNN,0.4,0.5,0.1
sample2,1,concat,concat,filter_DR_0.2,KNN,0.6,0.3,0.1
sample3,1,concat,concat,filter_DR_0.2,KNN,0.4,0.5,0.1
sample4,1,concat,concat,filter_DR_0.2,KNN,0.6,0.3,0.1
sample5,1,concat,filter_DR_0.2,KNN,concat,0.4,0.5,0.1
```

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
R                              4.2.3

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
BiocManager                    3.16
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
Yumei Li: ymli12@suda.edu.cn <br>
Junpeng Zhou: jpzhouzzz@stu.suda.edu.cn <br>
Keyao Zhu: kyzhu@stu.suda.edu.cn <br>

