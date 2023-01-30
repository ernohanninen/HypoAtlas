#Login:
ssh ernohan@kebnekaise.hpc2n.umu.se 

#Project folder:
cd /proj/nobackup/snic2022-22-293/ResearchProject/Training

#Dangpu
ssh bns631@dangpu01fl.unicph.domain

#https://sgn102.pages.ku.dk/a-not-long-tour-of-dangpu

#Move file from local to server:
scp file.py bns631@cortex:/home/bns631/

#Running jupyter on server
jupyter notebook --no-browser --port=8080

#Access jupyter
ssh -L 8080:localhost:8080 bns631@cortex


#Move dir from local to server
scp -r data ernohan@kebnekaise.hpc2n.umu.se:/proj/nobackup/snic2022-22-293/ResearchProject/Training

#Move dir from server to local
scp -r ernohan@kebnekaise.hpc2n.umu.se:/proj/nobackup/snic2022-22-293/ResearchProject/Training/Models/9_model  Models

scp -r ernohan@kebnekaise.hpc2n.umu.se:/proj/nobackup/snic2022-22-293/ResearchProject/Training/Models_2/2_model  Models_2



#Activate virtual environment from the directory environment located
source bin/activate

#Install python packages to env
pip install --no-cache-dir --no-build-isolation tensorflow

#Run pipeline
nextflow run main.nf -c /home/ernohanninen/SIMPLI/run/example_run.config 


