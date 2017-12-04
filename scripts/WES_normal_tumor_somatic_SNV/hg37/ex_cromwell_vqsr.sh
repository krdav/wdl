# load modules
module load moab torque
module load tools
module load oracle_jdk/1.8.0_144
module load gcc/6.2.0
module load R/3.2.5
module load mariadb/10.1.23

java -Dconfig.file=/home/projects/cu_10098/apps/cromwell.conf -jar /home/projects/cu_10098/apps/src/cromwell/target/scala-2.12/cromwell-30-c4b562f-SNAP.jar run WES_normal_tumor_somatic_SNV_vqsr.wdl -i WES_normal_tumor_somatic_SNV.inputs.json

