We use the `fork-and-pull` workflow for moose development.

## To contribute to the main moose code (moose-core):

   1. Fork the `development` branch in [moose-core](https://github.com/MooseNeuro/moose-core) repo on github
   2. Make the changes you want to contribute.
   3. Ensure that moose debug build works locally on your version (see build instructions to know how to build moose). Look at the compilation outputs and resolve all errors and warnings.
   5. Check that you are able to install it on your local system without errors.
   4. Check that all tests in the `tests` directory and its subdirectories pass.
   5. Submit your changes via a pull-request.

	
After this the repository admin will check your pull request and merge it to the main branch.

### How to add new classes to moose-core
	The core of moose is consists of a C++ framework using templates that is used for defining classes. You can find minimal examples of such class definitions in the `examples` subdirectory in the moose-core source code directory.
	
	
## How to contribute script/model examples	
	Please fork the [moose-examples](https://github.com/MooseNeuro/moose-examples) repository and add your example there. If it is a short standalone script, add it to the `snippets` subdirectory, if it is a complex model spanning multiple scripts/model files, please consider creating a new subdirectory with your example.

## Bug reports, feature suggestions are also contribution
	Please create a new issue on the [issue tracker](https://github.com/MooseNeuro/moose-core/issues) to report bugs or request features.
