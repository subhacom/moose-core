# Author: Subhasis Ray
# Date: Sun Jul  7

fmt_src = [
  'src/format.cc',
  'src/os.cc'
]

incdir = include_directories('include')

# This is to avoid deprecation warning from MSVC:
#   cpp_args: ['/D_SILENCE_ALL_MS_EXT_DEPRECATION_WARNINGS']
fmt_lib = static_library('fmt', fmt_src, include_directories: incdir, cpp_args: ['-D_SILENCE_ALL_MS_EXT_DEPRECATION_WARNINGS'])
