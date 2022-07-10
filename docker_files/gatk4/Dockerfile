# Set the base image to Ubuntu
FROM ubuntu:16.04

#Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
        liblzma-dev \
        libz-dev \
        libbz2-dev \
        libclang-dev \
        libcurl4-openssl-dev \
        python-pip \
        python2.7 \
        python2.7-dev \
        openjdk-8-jre && rm -rf /var/lib/apt/lists/* \
	&& apt-get clean

################## BEGIN INSTALLATION ######################
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 -O /tmp/htslib-1.15.1.tar.bz2\
    && tar xjf /tmp/htslib-1.15.1.tar.bz2 -C /tmp\
    && rm /tmp/htslib-1.15.1.tar.bz2\
    && cd /tmp/htslib-1.15.1\
    && make\
    && make install\
    && cd /\
    && rm -rf /tmp/htslib-1.15.1

RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip -O /tmp/gatk-4.2.6.1.zip\
	&& unzip /tmp/gatk-4.2.6.1.zip -d /opt/ \
	&& rm /tmp/gatk-4.2.6.1.zip -f \
	&& cd /opt/gatk-4.2.6.1 \
	&& ./gatk --list

ENV PATH="/opt/gatk-4.2.6.1/:${PATH}"
