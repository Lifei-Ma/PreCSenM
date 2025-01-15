# PreCSenM

**Cellular senescence** is characterized by a stable proliferation arrest triggered by various stressors. The **PreCSenM (Predictive Cellular Senescence Model)** enables the detection of senescence levels in samples from multiple cell types using gene expression data.



**Cellular Senescence Score Calculation**

**Description**

The PreCSenM function calculates Cellular Senescence (CS) scores from gene expression profiling data. It requires users to provide significant CS genes for analysis and supports multiple samples.



**Requirements and Input**

​	1.	**Gene Expression Matrix**

​	•	Upload a tab-separated value (TSV, TXT or CVS) file containing normalized gene expression data (e.g., TPM or FPKM values).

​	•	The file format should include:

​	•	The **first column**: Gene symbols.

​	•	Each **subsequent column**: Sample data with sample names in the header row and expression values below.

​	2.	**Normalization**

​	•	Ensure gene expression values are normalized (e.g., scaling).

​	•	The matrix should have **gene symbols as row names** and **sample IDs as column names**.

​	3.	**Significant Genes Data**

​	•	Users must load a predefined set of significant CS genes (sig.genes.Rdata) before running the function.



**Instructions for Use**

​	1.	Prepare your input file following the described format.

​	2.	Use the same quantification method for gene expression as the training data for accurate results.

​	3.	An example gene expression file is available for reference and download.



**Output**

The function returns a **dataframe** with the following columns:

​	•	**Sample_ID:** The IDs of the input samples.

​	•	**score.01:** Predicted senescence score:

​	•	**≥ 0.5 or median:** Classified as CS-high.

​	•	**< 0.5 or median:** Classified as CS-low.



**Example Usage**

\# Example usage of the PreCSenM function

\# Ensure to load normalized expression matrix and CS significant genes before running.

CS.score = PreCSenM(sig.genes = sig.genes, ExprMat = your_expression_matrix)



**Export Information**

​	•	**Export:** This function can be exported for use in R-based pipelines.

​	•	For any further assistance, you can review the provided example file. If you have any questions, feel free to reach out via email at lifei_ma@@126.com!