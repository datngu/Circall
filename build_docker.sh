## BUILD DOCKER AUTOMATICLY

docker build -t circall_dev:v0.0.1 -f Dockerfile .

# v1.0.0
docker tag circall_dev:v0.0.1 ndatth/circall:v1.0.1
docker push ndatth/circall:v1.0.1
echo DONE

## change bin data



## testing
cd /sigma4/data/genome_ref_GRCh37.75

# pull circall docker image:
docker pull ndatth/circall:v1.0.1

# test your docker
docker run --rm ndatth/circall:v1.0.1 Circall.sh

# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 createSqlite.R \
        data/Homo_sapiens.GRCh37.75.gtf \
        data/Homo_sapiens.GRCh37.75.sqlite

docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 buildBSJdb.R \
        gtfSqlite=data/Homo_sapiens.GRCh37.75.sqlite \
        genomeFastaFile=data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
        bsjDist=250 maxReadLen=150 \
        output=data/Homo_sapiens.GRCh37.75_BSJ_sequences.fa


docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh37.75.cdna.all.fa \
        -o data/IndexTranscriptome


docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh37.75_BSJ_sequences.fa \
        -o data/IndexBSJ



#wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_1.fasta.gz
#wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_2.fasta.gz

docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 Circall.sh \
    -genome data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    -gtfSqlite data/Homo_sapiens.GRCh37.75.sqlite \
    -txFasta data/Homo_sapiens.GRCh37.75.cdna.all.fa \
    -txIdx data/IndexTranscriptome \
    -bsjIdx data/IndexBSJ \
    -dep Circall/Data/Circall_depdata_human.RData \
    -read1 data/read1.fastq.gz \
    -read2 data/read2.fastq.gz \
    -p 20 \
    -tag testing_sample \
    -c FALSE \
    -td FASLE \
    -o data/real_data_sample12


##################################################
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_1.fasta.gz
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_2.fasta.gz

docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 Circall.sh \
    -genome data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    -gtfSqlite data/Homo_sapiens.GRCh37.75.sqlite \
    -txFasta data/Homo_sapiens.GRCh37.75.cdna.all.fa \
    -txIdx data/IndexTranscriptome \
    -bsjIdx data/IndexBSJ \
    -dep Circall/Data/Circall_depdata_human.RData \
    -read1 data/sample_01_1.fasta.gz \
    -read2 data/sample_01_2.fasta.gz \
    -p 12 \
    -tag testing_sample \
    -c FALSE \
    -td FASLE \
    -o data/Testing_out_v1.0.1




## testing hg38
cd /sigma4/data/genome_ref_hg38

wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz -O Homo_sapiens.GRCh38.106.gtf.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O Homo_sapiens.GRCh38.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


gunzip Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


# pull circall docker image:
docker pull ndatth/circall:v1.0.1

# test your docker
docker run --rm ndatth/circall:v1.0.1 Circall.sh


# assumming you are running unix and $PWD: is the path to directory that is mounted to docker containter as /data:
docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 createSqlite.R \
        data/Homo_sapiens.GRCh38.106.gtf \
        data/Homo_sapiens.GRCh38.106.sqlite


docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 buildBSJdb.R \
        gtfSqlite=data/Homo_sapiens.GRCh38.106.sqlite \
        genomeFastaFile=data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        bsjDist=250 maxReadLen=150 \
        output=data/Homo_sapiens.GRCh38.106_BSJ_sequences.fa


docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh38.cdna.all.fa \
        -o data/IndexTranscriptome


docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 TxIndexer \
        -t data/Homo_sapiens.GRCh38.106_BSJ_sequences.fa \
        -o data/IndexBSJ


wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_1.fasta.gz
wget https://github.com/datngu/Circall/releases/download/v0.1.0/sample_01_2.fasta.gz

docker run --rm -v $PWD:/data ndatth/circall:v1.0.1 Circall.sh \
    -genome data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -gtfSqlite data/Homo_sapiens.GRCh38.106.sqlite \
    -txFasta data/Homo_sapiens.GRCh38.cdna.all.fa \
    -txIdx data/IndexTranscriptome \
    -bsjIdx data/IndexBSJ \
    -dep Circall/Data/Circall_depdata_human.RData \
    -read1 data/sample_01_1.fasta.gz \
    -read2 data/sample_01_2.fasta.gz \
    -p 4 \
    -tag testing_sample \
    -c FALSE \
    -td FASLE \
    -o data/Testing_out














