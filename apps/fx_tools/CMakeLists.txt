cmake_minimum_required (VERSION 2.6)
project (sandbox_fx_tools_apps_fx_tools)

seqan_add_executable(fx_convert fx_convert.cpp)
seqan_add_executable(fx_faidx fx_faidx.cpp)
seqan_add_executable(fx_sak fx_sak.cpp)

seqan_add_executable(fx_read_fasta fx_read_fasta.cpp)
seqan_add_executable(fx_read_fastq fx_read_fastq.cpp)
seqan_add_executable(fx_write_fasta fx_write_fasta.cpp)
seqan_add_executable(fx_write_fastq fx_write_fastq.cpp)
seqan_add_executable(fx_fastq_stats fx_fastq_stats.cpp)

seqan_add_executable(fx_sam_coverage fx_sam_coverage.cpp)
