#load modules
module purge
module load plink2/2.00a3.7LM
module load gcc/12.3.0 python/3.12.2

#make plink files and .raw files

while read -r p

 do 
	#python3 updateRSID_rmMultAllelic.py "$p".vcf "$p"_newRS.vcf #make an rsID out of the row number remove multi-allelic sites

	#plink2 --vcf "$p"_newRS.vcf --make-bed --out "$p" #convert to vcf to plink	

	#plink2 --bfile "$p" -indep-pairwise 500kb 0.2 --out "$p" #linkage stats

	#shuf -n 10000 "$p".prune.in > "$p".prune_10k.in #get snps

	#plink2 --bfile "$p" --extract "$p".prune_10k.in --make-bed --out "$p"_ldPruned #ld prune

	#plink2 --bfile "$p"_ldPruned --export A --out "$p"_ldPruned #make GT input file

	echo "done $p" 

done < infiles.txt
