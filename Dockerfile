FROM python:3.9-alpine 

WORKDIR /opt/gcmpy

COPY requirements.txt  /opt/gcmpy/requirements.txt
COPY setup.py /opt/gcmpy/setup.py
COPY gcmpy /opt/gcmpy/gcmpy

RUN pip3 install -r /opt/gcmpyrequirements.txt

RUN python3 /opt/gcmpy/setup.py install
