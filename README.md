<h1> NGS-Analyzer Web App</h1>
<p> At this early stage of developement, <b><i>NGS-Analyzer Web Application</i></b> is capable of running RNA-seq pipline on bacteria samplse using BWA for indexing and mapping and FeatureCounts to produce count matrix and provide interactive output for downstream analysis using DESeq2.</P>
<hr>
<h3>File upload</h3>
<p>The app requires to provide amount of conditions (ex. Treated, Untreated) and submit samples for each conditions. It is important that all pair-ended files follow naming patterng Prefix follored by _r1.fastq or _r2.fastq (ex. Sample_r1.fastq, Sample_r2.fastq). The app also requires to provide sequenced genome and annotation files.</p>
<img src = 'https://drive.google.com/uc?export=view&id=1PQ4c23XhkKfMhXEXv3UDQ5hdazrfG1q6' alt = 'Visualization tools page: Scatter Plot'>

<hr>
<h3>Results panel</h3>
<p>After the pipline is done running it will take a user to results pannel for interactive analysis where user can explpre differentially expressed genes by adjust Log2FoldChange and p-value.</p>
<img src = 'https://drive.google.com/uc?export=view&id=1bey_rMbKHLpXHIRvSf0Ys-d3aF3HHZic' alt = 'Visualization tools page: Heatmap'>

