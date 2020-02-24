/* 
 * include requires tasks 
 */
include { index; quantification; fastqc; multiqc  } from './rnaseq-tasks.nf'
 
/* 
 * define the data analysis workflow 
 */
workflow rnaseqFlow {
    // required inputs
    take:
      transcriptome
      read_files
    // workflow implementation
    main:
      index(transcriptome)
      quantification(index.out, read_files)
      fastqc(read_files)
      multiqc( quantification.out.mix(fastqc.out).collect() )
}