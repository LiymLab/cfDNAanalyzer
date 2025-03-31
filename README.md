# cfDNAanalyzer
cfDNAanalyzer (<ins>c</ins>ell-<ins>f</ins>ree <ins>DNA</ins> sequencing data <ins>analyzer</ins>) is a toolkit for cfDNA genomic sequencing data analysis which includes three core modules: <br>
(1) Feature Extraction, which extracts multiple genomic and fragmentomic features at whole-genome or specific genomic-region levels; 
(2) Feature Processing and Selection, allowing for refinement and optimization of extracted features, suporting multiple processing and selection methods; 
(3) Machine Learning Model Building, supporting the development of both single-modality and multiple-modality predictive models for disease detection and classification.<br>
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
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.1  <ins>C</ins>opy <ins>N</ins>umber <ins>A</ins>lterations (CNA) ([<ins>Adalsteinsson *et al, Nat Commun*, 2017</ins>](https://www.nature.com/articles/s41467-017-00965-y))<br>
* Copy number alterations comprise deletions or amplifications of a particular region of the genome, with a size as low as a few kilobases up to entire chromosomes. 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2  <ins>E</ins>nd <ins>M</ins>otif frequency and diversity (EM) ([<ins>Lee *et al, PNAS*, 2018</ins>](https://www.pnas.org/doi/abs/10.1073/pnas.1815031116) ; [<ins>Jiang *et al, Cancer Discov*, 2020</ins>](https://doi.org/10.1158/2159-8290.CD-19-0622))
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
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.5  <ins>E</ins>nd <ins>M</ins>otif frequency and diversity for <ins>R</ins>egions (EMR)([<ins>Serpas *et al, PNAS*, 2018</ins>](https://www.pnas.org/doi/abs/10.1073/pnas.1815031116) ; [<ins>Jiang *et al, Cancer Discov*, 2020</ins>](https://doi.org/10.1158/2159-8290.CD-19-0622))
* We introduced end motif frequency and diversity for regions, which is defined as the frequency and diversity of all 256 4-mer end motifs for each region.<br>
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.6  <ins>F</ins>ragmentation <ins>P</ins>rofile for <ins>R</ins>egions (FPR) ([<ins>Cristiano *et al, Nature*, 2019</ins>](https://doi.org/10.1038/s41586-019-1272-6))
* We introduced fragmentation profile for regions, which is defined as the fraction of short cfDNA fragments (100–150 bp) to long cfDNA fragments (151–220 bp) for each region.

### 3.  Features for transcription start sites (TSSs): 
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.1 <ins>P</ins>romoter <ins>F</ins>ragmentation <ins>E</ins>ntropy (PFE) ([<ins>Esfahani *et al, Nat Biotechnol*, 2022</ins>](https://doi.org/10.1038/s41587-022-01222-4))<br>
* Promoter fragmentation entropy quantifies the diversity of cfDNA fragment lengths around the TSSs of genes. It is calculated by a modified Shannon index for lengths of cfDNA fragment where both ends fell within ±1 kb of the TSS. Then this cfDNA entropy measure is adjusted using a Dirichlet-multinomial model for normalization.
#### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3.2 <ins>TSS</ins> <ins>C</ins>overage (TSSC) ([<ins>Ulz *et al, Nat Genet*, 2016</ins>](https://www.nature.com/articles/ng.3648))
* TSS coverage refers to the cfDNA sequencing read coverage in the regions surrounding TSSs.

### Supported feature processing methods and machine learning models

### 1. Feature processing 

#### 1.1 Format transformation and missing value filter
* We initiate the preprocessing by loading configuration settings and sample labels to establish a structured approach. The configuration files are generated according to the specific features that users want to calculate.

* We first convert all input files for each sample into the corresponding format as described below. After formatting, we concatenate the formatted CSV files for each sample by merging them column-wise. Finally, we clean the dataset by removing any columns that contain missing values, ensuring the data is ready for downstream analysis. 

  ```shell
  sample,label,column1,column2,column3,...,columnN
  sample1,1,column11,column12,column13,...,column1N
  sample2,0,column21,column22,column23,...,column2N
  ...
  sampleX,0,columnX1,columnX2,columnX3,...,columnXN
  ```

 > The formatted CSV file consists of rows representing different samples, where each row represent one sample. The first column holds the sample identifiers (e.g., `sample1`, `sample2`, ..., `sampleX`), followed by a label column that indicates the sample's classification (e.g., `1` for positive, `0` for negative). 
 >
 > After the label, columns `column1` to `columnN` contain feature data for each sample. Each feature column  represents a specific feature ID corresponding to the sample in that row.

#### 1.2 Variance filter

* In the variance filtering step, we begin by assessing the features to identify those with variance **lower than a variance threshold (0.15 as default)**. By eliminating these low-variance features, our aim is to retain only those that have variability in our analysis. This step is **optional**, enabling users to tailor their processing steps to the specific needs of their dataset and analysis objectives.

#### 1.3 Z-score standardization

* We then perform **Z-score standardization** on each feature separately. This process transforms the data so that it has a mean of zero and a standard deviation of one. Standardizing the data in this way ensures that all features are on a comparable scale, which is important for many machine-learning techniques. Similar to the previous steps, this Z-score standardization is **optional**, providing users with the flexibility they need.

### 2. Feature selection

Feature selection is an important step in machine leanring. It involves removing irrelevant or redundant features based on importance rankings, which can help reduce model complexity and enhance the model's ability to make accurate predictions. Here we provide **four methods (filter, embedded, wrapper, and hybrid methods)** for feature selection ( [<ins>Pudjihartono *et al, Frontiers in Bioinformatics*, 2022</ins>](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9580915/)), as follows:

#### 2.1 Filter methods

Filter methods apply statistical tests to assess and rank features for their relevance to the target variable, independent of the classification algorithm. These methods can be executed by computing the correlation between features and the target variable, subsequently eliminating features with low scores.

| Methods                                           | Alias | Brief Introduction                                           | Source                                                      |
| ------------------------------------------------- | ----- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Missing Value Ratio                               | MR    | The Missing Value Ratio method removes features with a high percentage of missing data, ensuring that only features with sufficient information content are retained for further analysis. | [pandas](https://pandas.pydata.org/)                                                       |
| Information Gain                                  | IG    | Information Gain measures the reduction in uncertainty of the target variable when a specific feature is known. Higher information gain indicates a more informative feature for predicting the target variable. | [scikit-learn](https://scikit-learn.org/stable/index.html)                                                 |
| Chi-square Test                                   | CHI   | The Chi-square Test evaluates the independence between categorical features and the target variable. Features that show significant association with the target variable are selected based on their Chi-square statistics. | [scikit-learn](https://scikit-learn.org/stable/index.html)                                                 |
| Fisher's Score                                    | FS    | Fisher’s Score ranks features according to their ability to distinguish between classes by comparing the mean differences across classes relative to their variances. | [scikit-feature](https://github.com/jundongl/scikit-feature/tree/master) |
| FCBF Fast Correlation Filter                      | FCBF  | The FCBF method quickly identifies relevant features by analyzing their correlation with the target variable while removing redundant features that contribute little additional information. | [FCBF_module](https://github.com/SantiagoEG/FCBF_module)     |
| Permutation Importance                            | PI    | Permutation Importance assesses the significance of each feature by measuring the drop in model performance when the feature’s values are randomly permuted. Features causing significant degradation in performance are deemed important. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master)                                                  |
| Correlation Coefficient                           | CC    | Correlation Coefficient quantifies the linear relationship between features and the target variable. Features with higher correlation coefficients are considered more predictive and are selected for modeling. | [pandas](https://pandas.pydata.org/)                                                       |
| Low Variance Filter                               | LVF   | The Low Variance Filter eliminates features with minimal variability across samples, as such features are unlikely to be useful for distinguishing between classes. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master)                                                 |
| Mean Absolute Difference                          | MAD   | Mean Absolute Difference evaluates the average difference in feature values between classes. Features with larger differences are prioritized for greater discriminatory power. | [pandas](https://pandas.pydata.org/)                                                       |
| Dispersion Ratio                                  | DR    | Dispersion Ratio assesses the spread or variability of a feature across samples. Features with distinct dispersion patterns across classes are selected as they are likely to contribute to better classification. | [numpy](https://numpy.org/)                                                        |
| Mutual Information                                | MI    | Mutual Information measures the amount of information a feature shares with the target variable. Features with higher mutual information values are more relevant for model inclusion. | scikit-learn                                                 |
| ReliefF                                           | RLF   | ReliefF identifies features that consistently differentiate between instances of different classes by considering the nearest neighbors of each instance. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
| SURF                                              | SURF  | The SURF variant extends ReliefF by considering all instances within a predefined distance, enhancing its ability to capture more relevant features. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
| MultiSURF                                         | MSURF | MultiSURF dynamically adjusts the neighborhood size during selection, allowing it to capture complex, multi-feature interactions. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
| TuRF                                              | TRF   | TuRF iteratively removes the least informative features, refining the feature set to focus on those with the highest discriminative power. | [scikit-rebate](https://epistasislab.github.io/scikit-rebate/using/)                                                             |
#### 2.2 Embedded methods

Embedded methods incorporate feature selection directly into model training, allowing the learning algorithm to determine the most crucial features. Models such as Random Forest and LASSO regression effectively execute these methods, assigning importance to features throughout the learning phase.

| Methods                  | Alias      | Brief Introduction                                           | Source       |
| ------------------------ | ---------- | ------------------------------------------------------------ | ------------ |
| Lasso                    | LASSO      | Lasso performs feature selection by applying L1 regularization, which encourages sparsity in the model by driving the coefficients of less important features to zero. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| Ridge                    | RIDGE      | Ridge applies L2 regularization, which penalizes large coefficients, thereby reducing the impact of less important features while retaining all features in the model. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| ElasticNet               | ELASTICNET | ElasticNet combines the strengths of both Lasso (L1 regularization) and Ridge (L2 regularization), allowing for balanced feature selection with both sparsity and regularization. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| Random Forest Importance | RF         | Random Forest Importance ranks features based on their contribution to model accuracy. The method leverages the structure of decision trees to identify the most predictive features during the training process. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |

#### 2.3 Wrapper methods

Wrapper methods choose features according to their impact on a selected classifier’s performance, systematically seeking the optimal subset of features. This is done by evaluating various feature combinations and identifying the subset that enhances the model’s performance, typically using techniques such as recursive feature elimination.

| Methods                     | Alias | Brief Introduction                                                                                                                                                       | Source     |
|----------------------------|-------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------|
| Forward Feature Selection    | FFS   | Forward Feature Selection starts with an empty model, sequentially adding features that improve model performance until no further gain is achieved. This method builds the optimal feature set step by step. | [mlxtend](https://rasbt.github.io/mlxtend/)   |
| Backward Feature Selection    | BFS   | Backward Feature Selection begins with all features included and iteratively removes the least significant ones. The process continues until the removal of further features results in performance degradation. | [mlxtend](https://rasbt.github.io/mlxtend/)    |
| Exhaustive Feature Selection   | EFS   | Exhaustive Feature Selection evaluates all possible feature combinations, ensuring the identification of the best subset. Though computationally intensive, it provides the most thorough feature selection. | [mlxtend](https://rasbt.github.io/mlxtend/)    |
| Recursive Feature Elimination   | RFE   | Recursive Feature Selection recursively eliminates the least important features, refining the model by retaining only those that contribute the most to prediction accuracy. | [scikit-learn](https://github.com/jundongl/scikit-feature/tree/master) |
| Boruta Algorithm            | BOR   | The Boruta Algorithm iteratively compares the importance of real features with shuffled versions, selecting all features that consistently prove to be more informative than random noise. | [boruta](https://www.jstatsoft.org/article/view/v036i11)     |

#### 2.4 Hybrid methods

Hybrid approaches combine filter with wrapper or embedded methods to leverage their unique advantages. For instance, features are initially filtered using statistical tests, followed by a wrapper method that enhances the selection by examining interactions with the classifier.

| Methods                  | Alias                            | Brief Introduction                                           |
| ------------------------ | -------------------------------- | ------------------------------------------------------------ |
| Filter-Wrapper | FW | Filter methods are firstly used to reduce the feature set based on relevance, and then wrapper methods refine the selection to identify the optimal features for the model. |
| Filter-Embedded | FE | Filter methods are firstly used to reduce the feature set based on relevance, and then embedded methods refine the selection to identify the optimal features for the model. |

### 3. Single modality machine learning model

In this section, we will use supervised machine learning to analyze the **individual feature type** of all samples. The classifiers in the toolkit include **Random Forest, K-Nearest Neighbors, Gaussian Naive Bayes, Logistic Regression, Support Vector Machine, and eXtreme Gradient Boosting**. The goal of this approach is to assess how effective each classifier is in accurately identifying disease status based on individual feature type. Our methods are applicable to both **binary and multi-class classification** scenarios.

### 4. Multiple modality machine learning model

Then, we aim to use **multiple feature types** to improve the model performance. Our toolkit include **Concatenation-based, Model-based, and Transformation-based methods** ([Reel *et al, Biotechnology Advances*, 2021](https://pubmed.ncbi.nlm.nih.gov/33794304/)), applicable to **binary and multi-class classification** tasks, to maximize model effectiveness across different scenarios. 

#### 4.1 Concatenation-based methods

Concatenation-based methods combine features from different sources by directly appending them into a single feature vector. This straightforward approach allows models to access all available data but relies heavily on the subsequent learning algorithm to handle the potentially high-dimensional space effectively.

#### 4.2 Model-based methods

Model-based integration methods create multiple intermediate models for the different omics data and then build a final model from various intermediate models ([Reel *et al, Biotechnology Advances*, 2021](https://pubmed.ncbi.nlm.nih.gov/33794304/)). It boosts predictive accuracy by utilizing the strengths of individual models, enabling a deeper understanding of the interactions across different layers.

| Methods         | Brief Introduction                                           | Reference                                                    |
| --------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Average Voting  | Average Voting combines predictions from multiple models by calculating the mean of their outputs. This method can balance diverse feature impacts, reducing the effect of anomalies in individual features. |                                                              |
| Weighted Voting | Weighted Voting assigns weights to individual model predictions based on their performance, using either their proportion of the total AUC or inversely to their error rates. This strategy prioritizes more reliable features when integrating multiple feature sets, enhancing overall prediction accuracy. |                                                              |
| Majority Voting | Majority Voting selects the most common outcome among various model predictions. It ensures robustness by reducing the influence of sporadically erroneous model outputs. |                                                              |
| Stack Ensemble  | Stack Ensemble uses a meta-model to learn how to best combine the predictions of multiple base models. This method effectively integrates diverse feature-derived predictions to optimize final decision-making processes. | [<ins>Li *et al, Genome Med*, 2024</ins>](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01280-6) |

#### 4.3 Transformation-based methods

Transformation-based methods transform each of the omics datasets firstly into graphs or kernel matrices and then combine all of them into one before constructing a model ([Reel *et al, Biotechnology Advances*, 2021](https://pubmed.ncbi.nlm.nih.gov/33794304/)). This approach supports the integration of varied data by standardizing how they are represented.

| Category            | Methods                            | Brief Introduction                                           | Reference                                                    |
| ------------------- | ---------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Dimension reduction | PCA                                | PCA reduces dimensionality by transforming features into a set of linearly uncorrelated components. It simplifies multi-feature integration by focusing on components that explain the most variance, thereby enhancing model efficiency and clarity. | [Subramanian *et al, Bioinform Biol Insights*, 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7003173/) |
| Kernel matrix | Linear Kernel                      | The Linear Kernel measures direct linear relationships between features. It aids in multi-feature integration by emphasizing linear associations, making it suitable for linearly separable data. |                                                              |
| Kernel matrix | Polynomial Kernel                  | The Polynomial Kernel allows the capture of interactions between features at different degrees of complexity. It enhances multi-feature integration by providing a flexible framework to model nonlinear relationships in the data. |                                                              |
| Kernel matrix | Radial Basis Function (RBF) Kernel | The RBF Kernel maps features into a higher-dimensional space using a radial basis function, enabling effective classification of non-linearly separable datasets. This kernel is crucial for multi-feature integration as it can handle complex and non-linear interactions between features. |                                                              |
| Kernel matrix | Sigmoid Kernel                     | The Sigmoid Kernel projects data using a sigmoid function, similar to neural network activation functions. It supports multi-feature integration by transforming features into formats that highlight threshold-based classifications. |                                                              |
| Network | Similarity network fusion (SNF)    | SNF integrates multiple types of data by fusing similarity networks, reinforcing common structural features. It excels in multi-feature integration by constructing a holistic network view, revealing deep insights across combined datasets. | [Wang *et al, Nat Methods*, 2014](https://www.nature.com/articles/nmeth.2810) |

### Usage
```ruby
bash cfDNAanalyzer.sh -I <InputFile> -o <OutputDirectory> -F <Features> [Options]
```
#### Options: 
```
---------- Options for feature extraction ----------

----- General options -----
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
                          

---------- Options for downstream analysis ----------

----- Options for feature processing and selection -----
  --noFPS                               Skip the whole step of feature processing and selection. [optional]
  --varFrac                      FLOAT  Variance threshold as a proportion of the maximum variance. Features with variance below this threshold will be removed. Default: [0.15]
  --filterMethod                 STR    Filter methods employed for feature selection (MR IG CHI FS FCBF PI CC LVF MAD DR MI RLF SURF MSURF TRF).
                                        Methods should be set as a string separated by space, e.g., MR IG CHI. 
                                        Default: All available filter methods will be applied.
                                        The detailed description of each filter method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --filterFrac                   FLOAT  Fraction of features to retain when employing the filter method. Default: [0.2]
  --wrapperMethod                STR    Wrapper methods employed for feature selection (FFS BFS EFS RFS BOR).
                                        Methods should be set as a string separated by space, e.g., BOR RFS. 
                                        Default: All available wrapper methods will be applied.
                                        The detailed description of each wrapper method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --wrapperFrac                  FLOAT  Fraction of features to retain when employing the wrapper method. Default: [0.2]
  --embeddedMethod               STR    Embedded methods employed for feature selection (LASSO RIDGE ELASTICNET RF).
                                        Methods should be set as a string separated by space, e.g., LASSO RIDGE. 
                                        Default: All available embedded methods will be applied.
                                        The detailed description of each embedded method can be accessed at https://github.com/LiymLab/cfDNAanalyzer.
  --embeddedFrac                 FLOAT  Fraction of features to retain when employing the embedded method. Default: [0.2]
  --hybridType                   STR    Two methods used for hybrid method (FW/FE). Default: [FE]
  --hybridMethod1,hybridMethod2  STR    Subtype methods designated for method 1 and method 2 in the "--hybridType".
                                        Methods should be set as a string separated by space, e.g., BOR RFS.
                                        Default: All available subtype methods in method 1 and method 2 will be applied.
  --hybridFrac1,hybridFrac2       FLOAT Fraction of features to retain for method 1 and method 2 in the "--hybridType". Default: [0.2]
  
----- Options for machine learning model building -----
  --noMLMB                 Skip the whole step of machine learning model building. [optional]
  --classNum          INT  Number of classification categories (2 or more). Default: [2]
  --cvSingle          STR  Cross-validation method applied in single modality, options include leave-one-out (LOO) or K-fold (KFold). Default: [LOO]
  --nsplitSingle      INT  Number of folds designated for K-fold cross-validation in single modality. Default: [5]
  --classifierSingle  STR  Classifiers employed in single modality (KNN SVM RandomForest GaussianNB LogisticRegression XGB).
                           Classifiers should be set as a string separated by space, e.g., KNN XGB.
                           Default: All available classifiers will be applied.
  --cvMulti           STR  Cross-validation method applied in multiple modalities, options include leave-one-out (LOO) or K-fold (KFold). Default: [LOO]
  --nsplitMulti       INT  Number of folds designated for K-fold cross-validation in multiple modalities. Default: [5]
  --classifierMulti   STR  Classifiers employed in multiple modalities (KNN SVM RandomForest GaussianNB LogisticRegression XGB).
                           Classifiers should be set as a string separated by space, e.g., KNN XGB.
                           Default: All available classifiers will be applied.  
```

### Run the usage example
You can directly run cfDNAanalyzer by ```cfDNAanalyzer.sh``` provided in the ```cfDNAanalyzer/``` directory with the example files available at Zenodo (DOI:xxx).  
```ruby
bash cfDNAanalyzer.sh -I ./example/input/bam_input.txt -o ./example/output/ -F CNA,FP,OCF -g hg19 -b ./example/input/test.bed -f <reference.fa> --filterMethod MR IG CHI --wrapperMethod BOR --embeddedMethod LASSO RIDGE --classNum 2 --cvSingle LOO --classifierSingle KNN SVM RandomForest --cvMulti LOO --classifierMulti KNN SVM RandomForest > ./cfDNAanalyzer.log
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
│   │   └── std_filtered_[FeatureName].csv
│   └── Feature_Selection
│       ├── std_filtered_[FeatureName]_[embeddedMethod]_selectd.csv
│       ├── std_filtered_[FeatureName]_[filterMethod]_selectd.csv
│       └── std_filtered_[FeatureName]_[wrapperMethod]_selectd.csv
└── Machine_Learning
    ├── single_modality
    │   ├── single_modality_results.csv
    │   └── single_modality_roc_curves.pdf or [classifierSingle]_std_filtered_[featureName]_[filterMethod]/[wrapperMethod]/[embeddedFrac]_selected.csv_accuracy_bar_chart.pdf
    └── multiple_modalities
        ├── Concatenation_based
        │   ├── concatenation_based_classification_results.csv
        │   └── concatenation_based_roc_curves.pdf or [classifierMulti]_accuracy_bar_chart.pdf
        ├── Model_based
        │   ├── model_based_classification_results.csv
        │   └── model_based_roc_curves.pdf or [Model_based_Method]_accuracy_bar_chart.pdf
        └── Transformation_based
            ├── transformation_based_classification_results.csv
            └── transformation_based_classification_results.pdf or [classifierMulti]_[Transformation_based_Method]_accuracy_bar_chart.pdf
```
### Features
#### Copy Number Alterations (CNA)
```CNA.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain the estimated copy number for each bin. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each bin.
```r
sample,label,chr10_100000001_100100000,chr10_10000001_10100000,chr10_1000001_1100000,...
sample1,1,2,2,2,...
sample2,1,2,2,2,...
sample3,0,2,2,2,...
sample4,0,3,3,3,...
```

#### End Motif frequency and diversity of whole genome (EM)
```EM_motifs_frequency.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain motif frequency for input BAM file across the whole genome. Column name `[Motifs]` of these columns specifies different end motifs.
```r
sample,label,AAAA,AAAC,AAAG,...
sample1,1,0.018236894002843927,0.006713449465985273,0.008389349825354902,...
sample2,1,0.017000827317117054,0.006374204623569929,0.008194779752533857,...
sample3,0,0.01710512948689375,0.006104100121349162,0.008247394089259742,...
sample4,0,0.017105455019315003,0.0060242744221732985,0.008044839655444735,...
```
```EM_motifs_mds.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Column `MDS` after the label contain motif frequency diversity scores for input BAM file across the whole genome.
```r
sample,label,MDS
sample1,1,0.9656993783146016
sample2,1,0.9676102549001604
sample3,0,0.9699398580311792
sample4,0,0.9749328602806467
```

#### Fragmentation Profile of whole genome (FP)
```Fragmentation_Profile.txt``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain ratio of short to long fragments across the whole genome with a window of 100kb. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each window.
```r
sample,label,chr10_100000001_100100000,chr10_10000001_10100000,chr10_1000001_1100000,...
sample1,1,0.418811500602882,0.423061870382213,0.4643819178586,...
sample2,1,0.609854406338351,0.585338518996558,0.620208689660986,...
sample3,0,0.732437051581384,0.748931817389019,0.786052494499277,...
sample4,0,1.093014538658,1.17698788772542,1.25233575408901,...
```

#### Nucleosome Occupancy and Fuzziness (NOF)
```NOF_meanfuziness.csv``` and ```NOF_occupancy.csv``` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain average fuzziness value (```NOF_meanfuziness.csv```) or occupany value (```NOF_occupancy.csv```) for each region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,5.10545,15.2444,13.1614,...
sample2,1,14.3248,15.8161,21.7031,...
sample3,0,8.85557,15.4003,12.8956,...
sample4,0,91.5162,132.746,135.676,...
```
#### Nucleosome Profile (NP)
```NP_amplitude.csv```, ```NP_central_coverage.csv``` and ```NP_mean_coverage.csv``` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain nucleosome peak amplitude (```NP_amplitude.csv```), central coverage (```NP_central_coverage.csv```), and the mean coverage (```NP_mean_coverage.csv```) for each site set in the input directory. Column name `[site].10000.txt` of these columns specifies different sites set.
```r
sample,label,AHR.10000.txt,AR.10000.txt,ARID3A.10000.txt,...
sample1,1,0.98707,1.00713,0.97745,...
sample2,1,0.99656,1.01126,0.98364,...
sample3,0,0.97693,1.00172,0.97278,...
sample4,0,0.98448,1.01047,0.98699,...
```
#### Windowed Protection Score (WPS)
```WPS_long.csv``` and ```WPS_short.csv``` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain WPS for long fragments (```WPS_long.csv```), WPS for short fragments (```WPS_short.csv```) for each region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,-0.453047,-1.16434,-0.827173,...
sample2,1,-1.12138,-0.908591,-1.54246,...
sample3,0,-0.934066,-0.972028,-0.815684,...
sample4,0,-7.13536,-8.07143,-7.43057,...
```

#### Orientation-aware CfDNA Fragmentation (OCF)
```OCF.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain OCF values for each region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,0,434.782608695652,0,...
sample2,1,0,0,44.642857142857,...
sample3,0,0,0,0,...
sample4,0,0,-270.564042303173,52.9100529100527,...
```

#### End Motif frequency and diversity for regions (EMR)
Our toolkit outputs two types of EMR features at different levels:
```EMR_aggregated_motif_frequency.csv```, ```EMR_aggregated_mds.csv```, ```EMR_region_motif_frequency.csv``` and ```EMR_region_mds.csv``` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification.
&nbsp;<br>
(1) Motif frequency and diversity for all input regions aggregated together.
```EMR_aggregated_motif_frequency.csv```<br>
```r
sample,label,AAAA,AAAC,AAAG,...
sample1,1,0.013650548983279928,0.005704338714331273,0.007591231453780619,...
sample2,1,0.01251368286624472,0.005531035547384818,0.007236131187111442,...
sample3,0,0.012409680749096808,0.005262478802624788,0.007378376219617095,...
sample4,0,0.012241659095563946,0.005127588398094708,0.006955051282185939,...
```
```EMR_aggregated_mds.csv```<br>
```r
sample,label,MDS,...
sample1,1,0.9656993783146016,...
sample2,1,0.9676102549001604,...
sample3,0,0.9699398580311792,...
sample4,0,0.9749328602806467,...
```
&nbsp;<br>
(2) Motif frequency and diversity for each region separately.
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
```FPR_fragmentation_profile_regions.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain ratio of short to long fragments for each input region. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each region.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100190117_100192118,...
sample1,1,0.622159430128502,0.205120478226236,0.255837665883836,...
sample2,1,0.856960330876682,1.2464732274728,0.596067015830815,...
sample3,0,2.47597624814421,0.387795200405157,0.469562355131108,...
sample4,0,1.25467384738066,0.769168897029371,0.870592227235521,...
```


#### Promoter Fragmentation Entropy (PFE)
```PFE.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain PFE value for the 2kb surrounding TSS. Column name `[TSS]` of these columns specifies different TSS gene ID.
```r
sample,label,A1BG_1,A1CF_1,A2ML1_1,...
sample1,1,0.32419976285609,0.321489009591541,0.322693262670918,...
sample2,1,0.290434197788833,0.29136196112472,0.288107079740962,...
sample3,0,0.311475798027447,0.312773162689136,0.308581470806014,...
sample4,0,0.283555552446221,0.274671104972553,0.27558683027636,...
```

#### TSS Coverage (TSSC)
```TSSC_average_coverage.csv``` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain cfDNA fragment coverage value for each TSS surrounding regions. Column name `[chr]_[chrStart]_[chrEnd]` of these columns specifies the chromosome, start coordinate and end coordinate for each TSS surrounding regions.
```r
sample,label,chr10_100026951_100028952,chr10_100154064_100156065,chr10_100173939_100175940,...
sample1,1,0.4184664362135796,1.2892859019022593,0.9334820673085641,...
sample2,1,0.9917570255447303,1.0606034360844632,1.4275246769830265,...
sample3,0,0.49373285632500463,0.8862925235299335,0.7441468312554452,...
sample4,0,0.7481164091500803,1.01417140904693,1.0611052799677623,...
```

### Feature Processing and Selection

#### Feature Processing

`std_filtered_[FeatureName].csv` consists of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain **processed** feature data for each sample.
```r
Sample,label,chr16_27560234_27562235,chr5_150900709_150902710
sample1,1,-0.5046462083397414,1.3816213426919597
sample2,1,0.119404026538408,0.9470529008767274
sample3,0,-0.5046462083397414,-0.2443601515770356
```

#### Feature Selection

`std_filtered_[FeatureName]_[embeddedMethod]_selectd.csv`,`std_filtered_[FeatureName]_[filterMethod]_selectd.csv` and `std_filtered_[FeatureName]_[wrapperMethod]_selectd.csv` consist of rows representing different samples. The `sample` column holds the sample's file name, followed by a `label` column that indicates the sample's classification. Columns after the label contain **selected** feature data for each sample.

```r
Sample,label,chr16_27560234_27562235,chr5_150900709_150902710
sample1,1,-0.5046462083397414,1.3816213426919597
sample2,1,0.119404026538408,0.9470529008767274
sample3,0,-0.5046462083397414,-0.2443601515770356
```

### Two-class Machine Learning

#### Single Modality

`single_modality_results.csv` contains performance metrics for different classifiers and feature sets, detailing the classifier used, associated feature file, and various performance metrics such as accuracy, precision, recall, F1 score, AUC, average score, computation time, and memory usage (peak memory).

```r
Classifier,File,accuracy,precision,recall,f1,auc,average,Time_Taken,Memory_Usage_MB
KNN,file1.csv,0.49,0.47619047619047616,0.40816326530612246,0.43956043956043955,0.4489795918367347,0.4525787545787545,31.92484450340271,0.534846
SVM,file2.csv,0.59,0.58,0.5918367346938775,0.5858585858585857,0.5118047218887555,0.5719000084882437,0.7156319618225098,0.140312
SVM,file3.csv,0.24,0.23529411764705882,0.24489795918367346,0.24,0.765906362545018,0.34521968787515006,2.8905766010284424,0.683896
GaussianNB,file4.csv,0.42,0.4,0.3673469387755102,0.3829787234042553,0.3569427771108443,0.38545368785812195,0.7432112693786621,0.110663
```
`single_modality_roc_curves.pdf` compares the performance of different classifiers applied to various feature sets based on their AUC (Area Under the Curve) values. Each line represents the ROC curve for a specific classifier-feature combination.
![image](https://github.com/LiymLab/cfDNAanalyzer/blob/main/pics/two_single_modality_roc_curves.pdf)

#### Multiple Modalities
##### Concatenation based
`concatenation_based_classification_results.csv` contains performance metrics for different classifiers used in concatenation based method. The metrics include accuracy, precision, recall, F1 score, AUC, computation time, and memory usage (peak memory).

```r
Classifier,Accuracy,Precision,Recall,F1,AUC,Time_Taken,Memory_Usage_MB
KNN,0.4,0.3877551020408163,0.3877551020408163,0.3877551020408163,0.37775110044017607,0.07472038269042969,0.270047
SVM,0.53,0.5227272727272727,0.46938775510204084,0.49462365591397844,0.5480192076830731,0.10915851593017578,0.263
```

`concatenation_based_roc_curves.pdf` compares the performance of different classifiers used in concatenation based method based on their AUC (Area Under the Curve) values. Each line represents the ROC curve for a specific classifier.

<img src="/Users/zkey/Library/Application Support/typora-user-images/image-20240926161026506.png" alt="image-20240926161026506" style="zoom: 67%;" />

##### Model based
`model_based_classification_results.csv` contains performance metrics for different model based methods (every model based method can only use one classifier). The metrics include accuracy, precision, recall, F1 score, AUC, computation time, and memory usage (peak memory).

```r
Method,Accuracy,Precision,Recall,F1,AUC,Time_Taken,Memory_Usage_MB
majority,1.0,1.0,1.0,1.0,1.0,0.012170076370239258,0.0
average,1.0,1.0,1.0,1.0,1.0,0.010748863220214844,0.0
weighted_error,0.96,0.96,0.96,0.96,0.96,0.010942459106445312,0.0
stack_ensemble,1.0,1.0,0.96,0.9795918367346939,1.0,103.55750870704651,1.25
```

`model_based_roc_curves.pdf` compares the performance of different different model based methods based on their AUC (Area Under the Curve) values. Each line represents the ROC curve for a specific model based method.

<img src="/Users/zkey/Library/Application Support/typora-user-images/image-20240926161026506.png" alt="image-20240926161026506" style="zoom: 67%;" />

##### Transformation based
`transformation_based_classification_results.csv` contains performance metrics for different classifiers in different transformation based methods. The metrics include accuracy, precision, recall, F1 score, AUC, computation time, and memory usage (peak memory).

```r
Classifier,Transformation_method,Accuracy,Precision,Recall,F1_Score,AUC,Time_Taken,Memory_Usage_MB
KNN,PCA,0.94,1.0,0.88,0.9361702127659575,0.9935999999999999,8.09121036529541,2.723829
SVM,PCA,1.0,1.0,1.0,1.0,1.0,7.304638147354126,2.488549
RandomForest,PCA,1.0,1.0,1.0,1.0,1.0,33.66458058357239,2.243559
GaussianNB,PCA,0.94,0.9230769230769231,0.96,0.9411764705882353,0.9728,8.260015726089478,2.221695
LogisticRegression,PCA,1.0,1.0,1.0,1.0,1.0,7.151197910308838,2.590051
XGB,PCA,0.94,0.9230769230769231,0.96,0.9411764705882353,0.8864,88.6492531299591,2.043965
```

`transformation_based_classification_results.pdf` compares the performance of different classifiers in different transformation based methods based on their AUC (Area Under the Curve) values. Each line represents the ROC curve for a specific classifier in one transformation based method.

<img src="/Users/zkey/Library/Application Support/typora-user-images/image-20240926161026506.png" alt="image-20240926161026506" style="zoom: 67%;" />

### Multi-class Machine Learning

#### Single modality
`single_modality_results.csv` contains performance metrics for different classifiers and feature sets, detailing the classifier used, associated feature file, and various performance metrics such as accuracy, macro-precision, macro-recall, macro-f1, micro-precision, micro-recall, micro-f1, weighted auc, per class accuracy, computation time, and memory usage (peak memory).

```r
Classifier,File,accuracy,precision_macro,recall_macro,f1_macro,precision_micro,recall_micro,f1_micro,per_class_accuracy,Time_Taken,Memory_Usage_MB
RandomForest,file1.csv,0.23,0.1964573268921095,0.20023148148148145,0.19274220032840722,0.23,0.23,0.23,"{0: 0.25925925925925924, 1: 0.0, 2: 0.20833333333333334, 3: 0.3333333333333333}",31.92484450340271,0.534846
GaussianNB,file2.csv,0.23,0.2159722222222222,0.2122790404040404,0.21283431180691456,0.23,0.23,0.23,"{0: 0.2222222222222222, 1: 0.0625, 2: 0.2916666666666667, 3: 0.2727272727272727}",0.7156319618225098,0.140312
LogisticRegression,file3.csv,0.18,0.15425084175084175,0.15898569023569023,0.15506253006253007,0.18,0.18,0.18,"{0: 0.18518518518518517, 1: 0.0, 2: 0.20833333333333334, 3: 0.24242424242424243}",2.8905766010284424,0.683896
KNN,file4.csv,0.24,0.2187888198757764,0.22206439393939392,0.2186508695538028,0.24,0.24,0.24,"{0: 0.3333333333333333, 1: 0.0625, 2: 0.25, 3: 0.24242424242424243}",0.7432112693786621,0.110663
```

`[classifierSingle]_std_filtered_[featureName]_[filterMethod]/[wrapperMethod]/[embeddedFrac]_selected.csv_accuracy_bar_chart.pdf` compares the performance of different classifiers applied to various feature sets based on their overall accuracy and per-class accuracy. 

![image-20240926180340886](/Users/zkey/Library/Application Support/typora-user-images/image-20240926180340886.png)

#### Multiple modalities
##### Concatenation based
`concatenation_based_classification_results.csv` contains performance metrics for different classifiers used in concatenation based method. The metrics include accuracy, macro-precision, macro-recall, macro-f1, micro-precision, micro-recall, micro-f1, weighted auc, per-class accuracy, computation time, and memory usage (peak memory).

```r
Classifier,accuracy,precision_macro,recall_macro,f1_macro,precision_micro,recall_micro,f1_micro,per_class_accuracy,Time_Taken,Memory_Usage_MB
GaussianNB,0.6842105263157895,0.4603174603174603,0.5740740740740741,0.5103785103785103,0.6842105263157895,0.6842105263157895,0.6842105263157895,"{0: 0.8333333333333334, 1: 0.8888888888888888, 2: 0.0}",1.249962568283081,0.752197
KNN,0.8421052631578947,0.8917748917748917,0.75,0.7410256410256411,0.8421052631578947,0.8421052631578947,0.8421052631578947,"{0: 1.0, 1: 1.0, 2: 0.25}",2.4304635524749756,0.710494
SVM,0.8421052631578947,0.8917748917748917,0.75,0.7410256410256411,0.8421052631578947,0.8421052631578947,0.8421052631578947,"{0: 1.0, 1: 1.0, 2: 0.25}",2.6898982524871826,0.455536
LogisticRegression,0.9473684210526315,0.9523809523809524,0.9166666666666666,0.9267399267399267,0.9473684210526315,0.9473684210526315,0.9473684210526315,"{0: 1.0, 1: 1.0, 2: 0.75}",1.5078327655792236,0.0
RandomForest,0.8947368421052632,0.919047619047619,0.8333333333333334,0.8457040035987404,0.8947368421052632,0.8947368421052632,0.8947368421052632,"{0: 1.0, 1: 1.0, 2: 0.5}",10.915184259414673,0.0
XGB,0.47368421052631576,0.5,0.3796296296296296,0.3466666666666667,0.47368421052631576,0.47368421052631576,0.47368421052631576,"{0: 0.0, 1: 0.8888888888888888, 2: 0.25}",22.575932502746582,0.0
```

`[classifierMulti]_accuracy_bar_chart.pdf ` compares the performance of different classifiers used in concatenation based method based on their overall accuracy and per-class accuracy.

![image-20240926181529615](/Users/zkey/Library/Application Support/typora-user-images/image-20240926181529615.png)

##### Model based
`model_based_classification_results.csv` contains performance metrics for different model based methods (every model based method can only use one classifier). The metrics include accuracy, macro-precision, macro-recall, macro-f1, micro-precision, micro-recall, micro-f1, weighted auc, per-class accuracy, computation time, and memory usage (peak memory).

```r
method,accuracy,precision_macro,recall_macro,f1_macro,precision_micro,recall_micro,f1_micro,per_class_accuracy,Time_Taken,Memory_Usage_MB
majority,1.0,1.0,1.0,1.0,1.0,1.0,1.0,"{0: 1.0, 1: 1.0, 2: 1.0}",0.0015943050384521484,0.007029
average,0.9473684210526315,0.9666666666666667,0.9444444444444445,0.9521531100478469,0.9473684210526315,0.9473684210526315,0.9473684210526315,"{0: 0.8333333333333334, 1: 1.0, 2: 1.0}",0.005101680755615234,0.074313
weighted_error,1.0,1.0,1.0,1.0,1.0,1.0,1.0,"{0: 1.0, 1: 1.0, 2: 1.0}",0.005359172821044922,0.031187
stack_ensemble,0.7894736842105263,0.8974358974358975,0.6944444444444445,0.709090909090909,0.7894736842105263,0.7894736842105263,0.7894736842105263,"{0: 0.8333333333333334, 1: 1.0, 2: 0.25}",65.93589687347412,1.251697
```

`[Model_based_Method]_accuracy_bar_chart.pdf ` compares the performance of different model based methods based on their overall accuracy and per-class accuracy.

![image-20240926181529615](/Users/zkey/Library/Application Support/typora-user-images/image-20240926181529615.png)

##### Transformation based
`transformation_based_classification_results.csv` contains performance metrics for different classifiers in different transformation based methods. The metrics include accuracy, macro-precision, macro-recall, macro-f1, micro-precision, micro-recall, micro-f1, weighted auc, per-class accuracy, computation time, and memory usage (peak memory).

```r
Classifier,Transformation_method,accuracy,precision_macro,recall_macro,f1_macro,per_class_accuracy,Time_Taken,Memory_Usage_MB
KNN,PCA,0.9473684210526315,0.9666666666666667,0.9444444444444445,0.9521531100478469,"{1: 0.8333333333333334, 2: 1.0, 3: 1.0}",0.848419189453125,0.466756
SVM,PCA,0.9473684210526315,0.9666666666666667,0.9444444444444445,0.9521531100478469,"{1: 0.8333333333333334, 2: 1.0, 3: 1.0}",0.8906724452972412,0.393033,
RandomForest,PCA,0.9473684210526315,0.9666666666666667,0.9444444444444445,0.9521531100478469,"{1: 0.8333333333333334, 2: 1.0, 3: 1.0}",5.284086227416992,0.393033
GaussianNB,PCA,0.9473684210526315,0.9666666666666667,0.9444444444444445,0.9521531100478469,"{1: 0.8333333333333334, 2: 1.0, 3: 1.0}",0.8491766452789307,0.391665
LogisticRegression,PCA,0.9473684210526315,0.9666666666666667,0.9444444444444445,0.9521531100478469,"{1: 0.8333333333333334, 2: 1.0, 3: 1.0}",1.0212316513061523,0.391665
XGB,PCA,0.8947368421052632,0.9111111111111111,0.8611111111111112,0.879281537176274,"{0: 0.8333333333333334, 1: 1.0, 2: 0.75}",16.179036140441895,0.391665
```

`[classifierMulti]_[Transformation_based_Method]_accuracy_bar_chart.pdf ` compares the performance of different classifiers in different transformation based methods based on their overall accuracy and per-class accuracy.

![image-20240926181529615](/Users/zkey/Library/Application Support/typora-user-images/image-20240926181529615.png)

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

