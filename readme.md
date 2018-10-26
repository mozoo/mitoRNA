mitoRNA.py is a Python2.7 designed to assemble complete mitogenomes using RNA-Seq data and starting mitogenome(s) as reference. It uses an iterative approach, alternating reference mapping and de novo assembly. Additional tools required for running RNAtoMITO are Bowtie2 and Trinity. An help page can be racalled using the command -h.


Usage: MitoRNA.py [-h] -R FASTA -M MOD [-1 FASTQ1] [-2 FASTQ2] [-U SINGLE]
                  [-P THREADS] [-T TAXA] -O OUT [-C] [-N END] 

example: mitoRNA.py -I reference_mitogenomes.fasta -M paired -1 left.fastq -2 right.fastq -P 16 -T metazoa -O target_mitogenome

Required Arguments:
  -h, --help    show this help message and exit
  -R FASTA      Fasta file with the reference mitochondrial genome
  -M MOD        reads type: paired - unpaired
  -1 FASTQ1     mate 1 containing the raw reads
  -2 FASTQ2     mate 2 containing the raw reads
  -U SINGLE     Fastq containing the raw reads
  -O OUT        Output directory

Optional Arguments:
  -P THREADS    Number of threads, (1)
  -C            Full clean up (delete all intermediate files)
  -N END        maximum number of iterations
