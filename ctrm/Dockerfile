FROM ubuntu:devel

COPY . ./app
WORKDIR /app
RUN apt-get update && apt-get install -y apt-utils && apt-get install -y curl
RUN  apt-get update
RUN apt-get install nano
RUN   apt-get -y install build-essential && g++ --version
RUN   apt-get -y install libomp-dev
VOLUME /data
CMD bash
