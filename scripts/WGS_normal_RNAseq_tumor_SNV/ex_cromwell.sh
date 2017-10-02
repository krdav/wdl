
# -Dcall-caching.enabled=false
java -Dconfig.file=/home/projects/dp_00005/apps/cromwell.conf -jar /home/projects/dp_00005/apps/src/cromwell/target/scala-2.12/cromwell-30-c4b562f-SNAP.jar run WGS_normal_RNAseq_tumor_SNV.wdl -i WGS_normal_RNAseq_tumor_SNV.inputs.json

