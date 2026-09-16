# sad-docker

SAD docker container hosted on AWS ECR used to create containers in the cloud.

## Algorithm

As of version 2.0 this module runs [SADnm](https://github.com/kandread/SADnm) (Python,
numpy/scipy) instead of the original Julia `Sad.jl`. The output specification is
unchanged: `<reach_id>_sad.nc` with `valid`, `reach_id`, `A0`, `n`, `Qa`, `Q_u` and
`time_str`.

**`A0` and `n` are effective quantities, not surveyed or identified ones.** SADnm never
separately identifies Manning's `n` — roughness and slope are absorbed jointly into a
single level constant that is anchored to the monthly prior. Both are back-derived from
that constant, so bias in the prior climatology transfers directly into them, and `n`
additionally depends on the SWOT reach slope. Where a reach has no usable positive slope,
`n` is written as the fill value rather than being computed from a substituted one.

`Q_u` is the lognormal standard deviation implied by SADnm's predictive uncertainty.
Note that this uncertainty is currently known to be over-confident — a nominal 90%
interval captures roughly 70% empirically.

## publishing disclaimer

Please **do not share this** repository, algorithm code contained within, or containers with others. 

***Also note that all Docker containers are currently hosted in AWS. You may build the containers using the Dockerfiles located in each algorithm's directory.***

## Container

### Arguments

Each container takes the name of a reach file as an argument. This file contains a list of reaches or sets of reaches where each reach identifier or set is on one line.

### Input and output operations

#### FLPE
- EFS Input is mounted to the container at `/mnt/data/input`
- EFS Output is mounted to the container at `/mnt/data/output`

## AWS Environment variable

AWS_BATCH_JOB_ARRAY_INDEX is used to determine the reach identifier to process for all containers. The index value is used to retrieve a line in a file that lists all reach identifiers to be processed.

## deployment

There is a script to deploy the Docker container image and Terraform AWS infrastructure found in the `deploy` directory.

Script to deploy Terraform and Docker image AWS infrastructure

REQUIRES:

- jq (<https://jqlang.github.io/jq/>)
- docker (<https://docs.docker.com/desktop/>) > version Docker 1.5
- AWS CLI (<https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html>)
- Terraform (<https://developer.hashicorp.com/terraform/tutorials/aws-get-started/install-cli>)

Command line arguments:

[1] registry: Registry URI
[2] repository: Name of repository to create
[3] prefix: Prefix to use for AWS resources associated with environment deploying to
[4] s3_state_bucket: Name of the S3 bucket to store Terraform state in (no need for s3:// prefix)
[5] profile: Name of profile used to authenticate AWS CLI commands

Example usage: ``./deploy.sh "account-id.dkr.ecr.region.amazonaws.com" "container-image-name" "prefix-for-environment" "s3-state-bucket-name" "confluence-named-profile"`

Note: Run the script from the deploy directory.
