import subprocess,argparse,os,shutil,sys

#arguments parsing
parser=argparse.ArgumentParser(description="example: unicondria.py -I homo_sapiens.fasta -M paired -1 left.fastq -2 right.fastq -P 16 -T metazoa -O homosapiens_mito -C ")
#change the optional title
parser._optionals.title = "Required Arguments"
optional_args=parser.add_argument_group('Optional Arguments')
#input file
parser.add_argument('-R',dest='fasta',required=True,help='Fasta file with the reference mitochondrial genome')
parser.add_argument('-M',dest='mod',required=True,help='reads type: paired - unpaired')
parser.add_argument('-1',dest='fastq1',required=False,help='mate 1 containing the raw reads')
parser.add_argument('-2',dest='fastq2',required=False,help='mate 2 containing the raw reads')
parser.add_argument('-U',dest='single',required=False,help='Fastq containing the raw reads')
optional_args.add_argument('-P',dest='threads',required=False,help='Number of threads, (1)',default="1")
optional_args.add_argument('-T',dest='taxa',required=False,help='Taxa to filter the contigs (to avoid contaminations) (Metazoa)',default="Metazoa")
parser.add_argument('-O',dest='out',required=True,help='Output directory')
optional_args.add_argument('-C',dest='clean',required=False,help='Full clean up (delete all intermediate files)',action="store_true")
optional_args.add_argument('-N',dest='end',required=False,help='maximum number of iterations')
#parser.add_argument('-D',dest='nr',required=True,help='folder for nr database')
#optional_args.add_argument('-F',dest='threshold',required=False,help='threshold for contigs filtering (1)',default=1)
#parser.add_argument('-G',dest='gencode',required=True,help='Genetic code for blastx')


#def filter(script,fasta,gencode,taxon,threads,threshold,nr,iteration):		#function to filter trinity contigs
#	subprocess.call([script,taxon,fasta,threshold,threads,gencode,nr,iteration])	#it calls a bash script

args = parser.parse_args()

#here we get the path for the initial files, before changing the working directory

#script=os.path.dirname(os.path.realpath(sys.argv[0]))+"/filter.sh"

if args.fastq1:
	fastq1=os.path.abspath(args.fastq1)
	fastq2=os.path.abspath(args.fastq2)
elif args.single:
	single=os.path.abspath(args.single)

fasta=os.path.abspath(args.fasta)


if not os.path.exists(os.path.abspath(args.out)):
	os.makedirs(os.path.abspath(args.out))
	os.makedirs(os.path.abspath(args.out+"/temp"))
else:
	if not os.path.exists(os.path.abspath(args.out)+"/temp"):
		os.makedirs(os.path.abspath(args.out+"/temp"))
os.chdir(os.path.abspath(args.out)+"/temp")
#shutil.copyfile(script,"./filter.sh")
#os.chmod("./filter.sh",0o777)


check=False
#here we do the first database for bowtie2
subprocess.call(["bowtie2-build",fasta,"./database_0"])

a=0	#number of iterations
while check == False:
	filter_out="./"+str(args.taxa)+"_blastx_allmatching_"+str(a-1)+".fasta"
	#filter_out="Metazoa_1_blastx_atleastone_0.fasta"
	db="./database_"+str(a)
	reads="./aligned_"+str(a)
	alignment="./alignment_"+str(a)	#this will be used later to read the alignment perc.
	alignment_file=open("./alignment_"+str(a),"w")	#here we open the file to be writeble as stderr in the bowtie subprocess
	sam="./sam_"+str(a)
	tri_out="./trinity_"
	if a==0:	#in the first iteration we use the database from the reference fasta and we run the first bowtie on the starting read set
		if args.mod == 'paired':
			subprocess.call(["bowtie2","-x",db,"--al-conc",reads,"-p",args.threads,"--very-sensitive-local","-N","1","-L","10","-1",fastq1,"-2",fastq2,"-S",sam],stderr=alignment_file)
			subprocess.call(["Trinity","--seqType","fq","--max_memory","60G","--no_normalize_reads","--min_contig_length","150","--CPU",args.threads,"--full_cleanup","--left",reads+".1","--right",reads+".2","--output",tri_out+str(a)])
			#print "Filtering Trinity Contigs"
			#filter("./filter.sh",tri_out+str(a)+".Trinity.fasta",str(args.gencode),args.taxa,str(args.threads),str(args.threshold),nr,str(a))
			alignment_file.close()
			pass
		elif args.mod == 'unpaired':
			subprocess.call(["bowtie2","-x",db,"--al",reads,"-p",args.threads,"--very-sensitive-local","-N","1","-L","10","-U",single,"-S",sam],stderr=alignment_file)
			subprocess.call(["Trinity","--seqType","fq","--max_memory","60G","--no_normalize_reads","--min_contig_length","150","--CPU",args.threads,"--full_cleanup","--single",reads,"--output",tri_out+str(a)])
			#print "Filtering Trinity Contigs"
			#filter("./filter.sh",tri_out+str(a)+".Trinity.fasta",str(args.gencode),args.taxa,str(args.threads),str(args.threshold),nr,str(a))
		alignment_file.close()
		pass
	else:	#starting from the second iteration
		#for f in os.listdir("./"):
		#	print f
		#os.system("bowtie2-build"+str(filter_out)+db)	#non funzionava in nessun modo con subprocess e ho usato os.system
		subprocess.call(["bowtie2-build",os.path.abspath(tri_out+str(a-1)+".Trinity.fasta"),db])	#we use tri_out as input to the bowtie2 mapping
		#subprocess.call(["bowtie2-build",tri_out+str(a-1)+".Trinity.fasta",db])
		if args.mod == "paired":
			subprocess.call(["bowtie2","-x",db,"--al-conc",reads,"-p",args.threads,"--very-sensitive-local","-N","1","-L","10","-1",fastq1,"-2",fastq2,"-S",sam],stderr=alignment_file)
			subprocess.call(["Trinity","--seqType","fq","--max_memory","60G","--no_normalize_reads","--min_contig_length","150","--CPU",args.threads,"--full_cleanup","--left",reads+".1","--right",reads+".2","--output",tri_out+str(a)])
			#print "Filtering Trinity Contigs"
			#filter("./filter.sh",tri_out+str(a)+".Trinity.fasta",str(args.gencode),args.taxa,str(args.threads),str(args.threshold),nr,str(a))
		elif args.mod == "unpaired":
			subprocess.call(["bowtie2","-x",db,"--al",reads,"-p",args.threads,"--very-sensitive-local","-N","1","-L","10","-U",single,"-S",sam],stderr=alignment_file)
			subprocess.call(["Trinity","--seqType","fq","--max_memory","60G","--no_normalize_reads","--min_contig_length","150","--CPU",args.threads,"--full_cleanup","--single",reads,"--output",tri_out+str(a)])
			#print "Filtering Trinity Contigs"
			#filter("./filter.sh",tri_out+str(a)+".Trinity.fasta",str(args.gencode),args.taxa,str(args.threads),str(args.threshold),nr,str(a))
		alignment_file.close()
		A=open(alignment).readlines()
		B=open("./alignment_"+str(a-1)).readlines()
		if float(A[len(A)-1].split()[0].strip("%"))-float(B[len(B)-1].split()[0].strip("%")) < 0.01:	#if the new alignment percentage is almost identical to the previous one
			print "REACHED PLATEAU AT "+str(a+1)+" ITERATION"
			if args.clean == True:
				shutil.copyfile(tri_out+str(a)+".Trinity.fasta","../Trinity_final.fasta")
				#shutil.copyfile(filter_out,"../"+str(args.taxa)+"_"+str(args.threshold)+"_blastx_allmatching_"+str(a)+".fasta")
				os.chdir("../")
				shutil.rmtree("./temp")
			check = True
		elif a == args.end:
			print "REACHED "+str(a+1)+"NTH ITERATION"
			check=True 
	a=a+1
	
#metazoa_1_blastx_atleastone_0.fasta
