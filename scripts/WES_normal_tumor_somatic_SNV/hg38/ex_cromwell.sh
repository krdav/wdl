

# java -Dconfig.file=/home/projects/cu_10098/apps/cromwell.conf -Dcall-caching.enabled=false -jar /home/projects/cu_10098/apps/src/cromwell/target/scala-2.12/cromwell-30-c4b562f-SNAP.jar run WES_normal_tumor_somatic_SNV.wdl -i WES_normal_tumor_somatic_SNV.inputs.json
java -Dconfig.file=/home/projects/cu_10098/apps/cromwell.conf -jar /home/projects/cu_10098/apps/src/cromwell/target/scala-2.12/cromwell-30-c4b562f-SNAP.jar run WES_normal_tumor_somatic_SNV.wdl -i WES_normal_tumor_somatic_SNV.inputs.json

