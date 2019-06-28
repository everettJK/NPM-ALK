#!/bin/sh
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_25.fq.gz -o ~/jp201805/salmon/WT_D6_4
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_26.fq.gz -o ~/jp201805/salmon/KD_D6_4
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_27.fq.gz -o ~/jp201805/salmon/TrpM_D6_4
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_28.fq.gz -o ~/jp201805/salmon/KD_D9_4
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_29.fq.gz -o ~/jp201805/salmon/TrpM_D9_4
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_30.fq.gz -o ~/jp201805/salmon/WT_D6_5
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_32.fq.gz -o ~/jp201805/salmon/KD_D6_5
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_33.fq.gz -o ~/jp201805/salmon/TrpM_D6_5
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_34.fq.gz -o ~/jp201805/salmon/WT_D9_5
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_35.fq.gz -o ~/jp201805/salmon/KD_D9_5
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_36.fq.gz -o ~/jp201805/salmon/TrpM_D9_5
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_37.fq.gz -o ~/jp201805/salmon/WT_D6_6
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_38.fq.gz -o ~/jp201805/salmon/KD_D6_6
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L01_39.fq.gz -o ~/jp201805/salmon/TrpM_D6_6
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L02_40.fq.gz -o ~/jp201805/salmon/WT_D9_6
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L02_41.fq.gz -o ~/jp201805/salmon/KD_D9_6
bsub -q fs02 -n 8 -M 24000 salmon quant -p 8 -i ~/GRCh38.p12/gencode28/salmonIndex/ -l A -r ~/jp201805/fastq/CL100021371_L02_42.fq.gz -o ~/jp201805/salmon/TrpM_D9_6

