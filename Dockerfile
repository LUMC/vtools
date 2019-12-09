FROM debian:buster-slim

RUN apt update && apt install python3-pip python3-dev cython3 zlib1g-dev python3-numpy libssl-dev python3-setuptools -y
RUN apt install -y libcurl4-openssl-dev
RUN apt install -y libbz2-dev liblzma-dev
COPY . vtools

RUN cd vtools && pip3 install .
