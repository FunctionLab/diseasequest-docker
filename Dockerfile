FROM centos:7

RUN yum install perl -y

ENV DQ_HOME=/dq/

RUN yum update -y

RUN mkdir -p $DQ_HOME/data $DQ_HOME/src $DQ_HOME/bin $DQ_HOME/outputs/networks

COPY ./bin/ $DQ_HOME/bin
COPY ./src/ $DQ_HOME/src


WORKDIR $DQ_HOME

COPY ./docker-entrypoint.sh /
ENTRYPOINT ["/docker-entrypoint.sh"]
