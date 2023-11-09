Assume we're starting from a clean repo and a fully set up conda environment::
        
    ```bash 
    ./bootstrap-conda
    mamba env create --file src/environment-dev-3.11.yml --name sage-dev
    conda activate sage-dev
    ```

(Note, that in the codespace environment you first have to delete the
already compiled files, e.g. with ``shopt -s globstar`` followed by ``rm src/**/*.so``
or ``for f in src/**/*.so ; do mv "$f" "$f.old"; done``
)

To compile and install the project in editable install, just use
    
    ```bash
    pip install --no-build-isolation --config-settings=builddir=builddir --editable .
    ```

Under the hood, pip invokes meson to configure and build the project.
We can also use meson directly as follows.

Now to configure the project, we need to run the following commands::

    ```bash
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
