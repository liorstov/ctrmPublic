# syntax=docker/dockerfile:1
FROM python:3
WORKDIR /
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
COPY rays_model.py rays_model.py
COPY start.sh start.sh
ENTRYPOINT /start.sh
