FROM python:3.8

COPY . .

WORKDIR /src
COPY requirements.txt .

RUN apt-get update 
RUN pip install -r requirements.txt

