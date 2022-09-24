FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:dd8f-main

RUN apt-get update --yes &&\
    apt-get install --yes --no-install-recommends wget software-properties-common dirmngr


RUN apt-get install -y curl vim default-jre-headless zlib1g zlib1g-dev unzip cmake
RUN python3 -m pip install cutadapt
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz &&\
    tar xvzf trim_galore.tar.gz &&\
    mv TrimGalore-0.6.6/trim_galore /bin &&\
    rm -rf TrimGalore-0.6.6

RUN curl -s https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O &&\
    unzip *.zip &&\
    chmod u+x ./FastQC/fastqc

RUN apt-get install -y liblzma-dev libncurses-dev libbz2-dev libssl-dev libcurl4-openssl-dev
RUN curl -L https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 -o samtools-1.16.tar.bz2 &&\
    tar -vxjf samtools-1.16.tar.bz2 &&\
    cd samtools-1.16 &&\
    make &&\
    make install

RUN curl -L https://tukaani.org/xz/xz-5.2.6.tar.gz -o xz-5.2.6.tar.gz &&\
    tar -xzvf xz-5.2.6.tar.gz &&\
    cd xz-5.2.6 &&\
    ./configure --enable-shared &&\
    make &&\
    make install &&\
    ldconfig

RUN wget https://sourceforge.net/projects/libpng/files/zlib/1.2.9/zlib-1.2.9.tar.gz/download &&\
    tar -xzvf download &&\
    cd zlib-1.2.9 &&\
    ./configure && make && make install

RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.5/bowtie2-2.4.5-linux-x86_64.zip &&\
    unzip bowtie2-2.4.5-linux-x86_64.zip &&\
    mv bowtie2-2.4.5-linux-x86_64 /bin

ENV PATH="${PATH}:/bin/bowtie2-2.4.5-linux-x86_64"

RUN python3 -m pip install --upgrade multiqc matplotlib numpy scipy lgenome cykhash macs3

RUN python3 -m pip install latch --upgrade
COPY wf/ /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
