samples = snakemake.params["samples"]
genome = snakemake.params["genome"]

f = open("./metadata/majiq_config_file.txt","w")
f.write("[info]\n")
f.write("readlen=150\n")
f.write("bamdirs=./files/bam_files\n")
f.write("genome="+genome+"\n")
f.write("strandness=None\n")
f.write("[experiments]\n")

for rep in samples:
	f.write(rep+'='+rep+'_Aligned.sortedByCoord.out\n')
