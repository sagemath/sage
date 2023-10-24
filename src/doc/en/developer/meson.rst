Assume we're starting from a clean repo and a fully set up conda environment::
        
    ```bash 
    ./bootstrap-conda
    mamba env create --file src/environment-dev-3.11.yml --name sage-dev
    conda activate sage-dev
    ```

Now to configure the project, we need to run the following commands::

    ```bash
    ./bootstrap
    meson setup builddir
    ```

This will create a build directory ``builddir`` that will hold the build artifacts.
To compile the project, run the following command::

    ```bash
    meson compile -C builddir
    ```

Installing the project is done with the following command::

    ```bash
    meson install -C builddir
    ```
