<h1> NGS-Analyzer Web App</h1>
<p> At this early stage of development, the <b><i>NGS-Analyzer Web Application</i></b> is capable of running an RNA-seq pipeline on bacterial samples, utilizing BWA for indexing and mapping, and FeatureCounts to generate a count matrix. The application also provides interactive output for downstream analysis with DESeq2.</p>
<hr>
<h3>File upload</h3>
<p>The app requires specifying the number of conditions (e.g., Treated, Untreated) and submitting samples for each condition. It is important that all paired-end files follow a naming pattern: a prefix followed by _r1.fastq or _r2.fastq (e.g., Sample_r1.fastq, Sample_r2.fastq). Additionally, the app requires providing sequenced genome and annotation files.</p>

<div align="center">

  ![ezgif com-video-to-gif-converter (2)](https://github.com/user-attachments/assets/4b1284ee-8150-48be-bec7-ff0947275aa5)

</div>

<hr>
<h3>Results panel</h3>
<p>Once the pipeline has finished running, the user will be directed to the results panel for interactive analysis. Here, users can explore differentially expressed genes by adjusting the Log2FoldChange and p-value.</p>

<div align="center">

  ![RESULTS-ezgif com-video-to-gif-converter](https://github.com/user-attachments/assets/b8362103-9732-415b-92e7-be4c4de0b755)


</div>
