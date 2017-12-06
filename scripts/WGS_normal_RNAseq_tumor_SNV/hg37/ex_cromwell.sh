# load modules
module load moab torque
module load tools
module load oracle_jdk/1.8.0_144
module load gcc/6.2.0
module load R/3.2.5
module load mariadb/10.1.23
module load star/2.5.3a

# -Dcall-caching.enabled=false
java -Dconfig.file=/home/projects/dp_00005/apps/cromwell.conf -jar /home/projects/dp_00005/apps/src/cromwell/target/scala-2.12/cromwell-30-c4b562f-SNAP.jar run WGS_normal_RNAseq_tumor_SNV.wdl -i WGS_normal_RNAseq_tumor_SNV.inputs.json

