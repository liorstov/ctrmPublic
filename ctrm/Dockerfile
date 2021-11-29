
FROM ubuntu:devel

COPY /build ./app/build
#COPY /testData ./app/testData
COPY ./script.sh ./app
WORKDIR /app
RUN  apt-get update &&  apt-get install -y libgomp1  && apt-get install python3.8
RUN apt-get install nano
RUN pip install numpy, pandas
CMD bash
