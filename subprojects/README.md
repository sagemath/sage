# Meson Subprojects

This directory contains subprojects used by the Meson build system. Subprojects are a way to include external dependencies in the Meson project and build them alongside the main Sagemath project. 

The metadata for each subproject is defined in a `*.wrap` file, which specifies the source files, version, and other relevant information about the subproject. Moreover, most subprojects have a corresponding folder under `packagefiles` that contains the necessary files for building and relevant patches.

For more detailed information, refer to the Meson documentation on [subprojects](https://mesonbuild.com/Subprojects.html) and [wrap files](https://mesonbuild.com/Wrap-dependency-system-manual.html).
