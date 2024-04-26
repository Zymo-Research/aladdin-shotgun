FROM nfcore/base:2.1

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/sourmash/bin:$PATH