#!/bin/bash
# To run
#               ./bamMove.sh All.DA_LowPass.Samples_Movelist	   #
#

while read line1
do

######## Variables to change #########

SRR=$line1

# Path to folder to store processed bams and bais
outputFolder='/path/to/mappingfolder/BamBaiOut'
# Path to output folder
baseFolder='/path/to/mappingfolder/'

echo "#!/bin/bash -l" > bamMove.sc;
echo "#SBATCH -A $proj_id" >> bamMove.sc;
echo "#SBATCH -p core -n 1" >> bamMove.sc;
echo "#SBATCH -J $SRR-bamMove" >> bamMove.sc;
echo "#SBATCH -t 10:00" >> bamMove.sc;

# load some modules
# echo "module load bioinfo-tools" >> bamMove.sc;
# echo "module load bcftools" >> bamMove.sc;

echo "SRR=$SRR" >> bamMove.sc;
echo "outputFolder=$outputFolder" >> bamMove.sc;
echo "baseFolder=$baseFolder" >> bamMove.sc;

# Move the bams and bais from per sample folders in one common folder
echo "cd $baseFolder/$SRR" >> bamMove.sc;
echo "mv $SRR.sorted.merged.bam $outputFolder" >> bamMove.sc;
echo "mv $SRR.sorted.merged.bam.bai $outputFolder" >> bamMove.sc;

sbatch bamMove.sc;

done<All.DA_LowPass.Samples_Movelist

#All.DA_LowPass.Samples_Movelist is just a list of samplesIDs.