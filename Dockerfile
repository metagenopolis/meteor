FROM continuumio/miniconda3

# set work dir
WORKDIR /app

# install conda packages
RUN conda create --name meteor -c conda-forge -c bioconda  -c aghozlane meteor && conda clean -ya

# add meteor environment to path
RUN echo "source activate meteor" > ~/.bashrc
ENV PATH /opt/conda/envs/meteor/bin:$PATH

ENTRYPOINT ["meteor"]

CMD [ "meteor -h" ]