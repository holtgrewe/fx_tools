cmake_minimum_required (VERSION 2.6)
project (sandbox_fx_tools_apps_fx_tools)

seqan_add_executable(fx_convert fx_convert.cpp)
seqan_add_executable(fx_faidx fx_faidx.cpp)
seqan_add_executable(fx_sak fx_sak.cpp)
seqan_add_executable(fx_renamer fx_renamer.cpp)

seqan_add_executable(fx_sam_coverage fx_sam_coverage.cpp)
