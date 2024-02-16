# Use an official Python runtime as the base image
FROM python:3.9
WORKDIR /usr/src/app
RUN useradd -m -s /bin/bash appuser
RUN chown appuser:appuser /usr/src/app
# Set the working directory in the container

COPY requirements.txt /usr/src/app/

RUN pip install --no-cache-dir -r requirements.txt
RUN jupyter nbextension enable  --py --sys-prefix graph_notebook.widgets
RUN python -m graph_notebook.static_resources.install
RUN python -m graph_notebook.nbextensions.install
RUN python -m graph_notebook.ipython_profile.configure_ipython_profile

USER appuser


EXPOSE 8888


CMD ["jupyter", "notebook", "--NotebookApp.kernel_manager_class=notebook.services.kernels.kernelmanager.AsyncMappingKernelManager", "--ip", "0.0.0.0", "/usr/src/app/"]