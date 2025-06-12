Steps to View Datasets

1. Load the docker image with:
   `docker load -i handle_dataset_views.tar.gz`
2. Put the dataset view files (*.sqlite) files in a directory, "views"
3. Make a directory "outputs"
4. Run the docker image to spawn the jupyter notebook.
    $ mkdir views; mv *sqlite views/
    $ mkdir outputs
    $ docker run -p 8888:8888 -v ./views:/workspace/views -v ./outputs:/workspace/outputs handle_dataset_views
    The `-p` flag exposed the port `8888` inside the docker image to the port by the same name externally. 
    The `-v` flag exposes a directory (in this case `./views`, so put your dataset views there) to a directory inside the docker image so that the jupyter notebook and access them.
    The ./outputs directory provides another shared directory that can be useful to pass output files.
    If using a M* MAC, the flag `--platform=linux/amd64` could be required.
5. Entering the URL that starts with `http://127.0.0.1:8888...` in a internet browser should lead to a jupyterlab instance.
6. Use the jupyter notebook to access the data. Because the notebook is running in the docker image, you'll need to put files that you want externally in the "output" directory created in step 3.