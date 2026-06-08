# HW3 Nextflow pipeline
## Profiles
- local:  for execution on the current machine and uses a local conda environment path passed with --conda_env.  
- cluster: for execution on remote cluster.  
- container: for execution of all processes in the Docker image from --container.

## Docker image

```
docker login 
docker build  -t mashakrasha/hw3-nextflow:latest .
docker run --rm  mashakrasha/hw3-nextflow:latest
```

```
docker pull mashakrasha/hw3-nextflow:latest
```

## Run examples
1) local
```
nextflow run main.nf \
    --reference ref_HV.fna
    -profile local \
    --samplesheet samples.csv \
    -stub-run
```
2) cluster
```
nextflow run main.nf \
  --reference ref_HV.fna
  -profile cluster \
  --samplesheet samples.csv 

```
3) container 
```
nextflow run main.nf \
  --reference ref_HV.fna
  -profile container \
  --samplesheet samples.csv 
``` 
## Example data 
You can find [here](https://disk.yandex.ru/d/GP9P5OE6j-elUQ)
