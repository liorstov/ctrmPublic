FROM ubuntu:devel

COPY /build ./app/build
#COPY /testData ./app/testData
COPY ./script.sh ./app
WORKDIR /app
RUN  apt-get update &&  apt-get install -y libgomp1
RUN apt-get install nano

CMD bash
