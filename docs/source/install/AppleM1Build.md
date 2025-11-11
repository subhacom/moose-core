# Building MOOSE on MacOS with Apple M1 CPU
- Install homebrew: https://brew.sh/
- Set up required development environment
  - Install command line tools for XCode
  - Install build dependencies by running these commands in a terminal
  ```
          brew install gsl
  ```

- Install anaconda/miniconda/micromamba/miniforge. For example, for micromamba, run
```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

  in command line
- Update micromamba: `micromamba self-update`
- Restart terminal and create an environment with necessary packages:

```
micromamba create -n moose numpy matplotlib vpython lxml meson ninja meson-python gsl setuptools pybind11[global] pkg-config -c conda-forge
```

- Activate the moose environment: `micromamba activate moose`

- Build and install using meson

```
cd moose-core
meson setup --wipe _build -Duse_mpi=false --buildtype=release
meson compile -vC _build
meson install -C _build
```

This will install moose at the standard location of the environment. To install at a specific location, you have to pass the full path with the `--prefix` argument during meson setup. For example, `_build_install` subdirectory of the current directory can be specified with:

```
meson setup --wipe _build --prefix=`pwd`/_build_install -Duse_mpi=false --buildtype=release
```

  - **Buildtype**
	If you want a developement build with debug enabled, pass `-Dbuildtype=debug` in the `meson setup`.


	```
	meson setup --wipe _build --prefix=`pwd`/_build_install -Duse_mpi=false -Dbuildtype=debug
	```

	You can either use `buildtype` option alone or use the two options `debug` and `optimization` for finer grained control over the build. According to `meson` documentation `-Dbuildtype=debug` will create a debug build with optimization level 0 (i.e., no optimization, passing `-O0 -g` to GCC), `-Dbuildtype=debugoptimized`  will create a debug build with optimization level 2 (equivalent to `-Ddebug=true -Doptimization=2`), `-Dbuildtype=release` will create a release build with optimization level 3 (equivalent to `-Ddebug=false -Doptimization=3`), and `-Dbuildtype=minsize` will create a release build with space optimization (passing `-Os` to GCC).

  - **Optimization level**

	To set optimization level, pass `-Doptimization=level`, where level can be `plain`, `0`, `g`, `1`, `2`, `3`, `s`.
