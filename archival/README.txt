Steps to View Datasets

Follow the steps below to view these datasets. Ensure Docker is installed and running before starting.

1. Load the Docker image with:
   `docker load -i handle_dataset_views.tar.gz`
2. Put the dataset view files (*.sqlite) files in a directory, "views"
    $ mkdir views; mv *sqlite views/
3. Make a directory "outputs"
    $ mkdir outputs
4. Run the Docker image to spawn the Jupyter notebook.
    $ docker run -p 8888:8888 -v ./views:/workspace/views -v ./outputs:/workspace/outputs handle_dataset_views
   The `-p` flag exposes the port `8888` inside the Docker image to the port by the same name externally. Change this number if you wish to use a different port.
   The `-v` flag exposes a directory (in this case `./views`, so put your dataset views there) to a directory inside the Docker image so that the Jupyter notebook and access them.
   The ./outputs directory provides another shared directory that can be useful for saving output files.
   If using a M* Mac, the flag `--platform=linux/amd64` could be required.
5. Entering the URL that starts with `http://127.0.0.1:8888...` in a internet browser should lead to a JupyterLab instance.
6. Use the Jupyter notebook to access the data. Because the notebook is running in the Docker image, you'll need to put files that you want externally in the "output" directory created in step 3.