// ==========================================================================
//                               FX Tools
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Renaming of FASTA/Q read identifiers.
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// Enum RenameSchema
// --------------------------------------------------------------------------

// Represent the renaming schema.

enum RenameSchema
{
    RENAME_SEQ,     // Rename to sequence.
    RENAME_NUMERIC  // Rename to counter.
};

// --------------------------------------------------------------------------
// Class FxRenamerOptions
// --------------------------------------------------------------------------

struct FxRenamerOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;
    
    // Path to FASTA/FASTQ file.
    seqan::CharString inPath;

    // Path to output file.
    seqan::CharString outPath;

    // The schema to use for renaming.
    RenameSchema renameSchema;

    FxRenamerOptions() : verbosity(1), renameSchema(RENAME_SEQ)
    {}
};

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxRenamerOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_renamer");
    setShortDescription(parser, "Renaming of FASTA and FASTQ files.");
    setVersion(parser, "0.1");
    setDate(parser, "August 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] -i \\fIOUT.fasta\\fP \\fB-o\\fP \\fIOUT.fasta\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] -i \\fIOUT.fastq\\fP \\fB-o\\fP \\fIOUT.fastq\\fP");
    addDescription(parser, "Rename sequences in FASTA and FASTQ files based on their sequence or an increasing counter.");

    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "I/O Options");
    addOption(parser, seqan::ArgParseOption("i", "in-path", "Path to input file.", seqan::ArgParseArgument::INPUTFILE, "FILE"));
    setRequired(parser, "in-path", true);
    setValidValues(parser, "in-path", "fasta fa fastq fq");
    addOption(parser, seqan::ArgParseOption("o", "out-path", "Path to the resulting file.  If omitted, result is printed to stdout.", seqan::ArgParseArgument::OUTPUTFILE, "FILE"));
    setRequired(parser, "out-path", true);
    setValidValues(parser, "out-path", "fasta fa fastq fq");

    addSection(parser, "Rename Options");
    addOption(parser, seqan::ArgParseOption("s", "rename-schema", "Select the renaming schema.", seqan::ArgParseArgument::STRING, "SCHEMA"));
    setValidValues(parser, "rename-schema", "sequence numeric");
    setDefaultValue(parser, "rename-schema", "sequence");
    
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getOptionValue(options.inPath, parser, "in-path");
        getOptionValue(options.outPath, parser, "out-path");

        seqan::CharString renameSchema;
        getOptionValue(renameSchema, parser, "rename-schema");
        if (renameSchema == "sequence")
            options.renameSchema = RENAME_SEQ;
        else
            options.renameSchema = RENAME_NUMERIC;

        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;
    }

    return res;
}

// ---------------------------------------------------------------------------
// Function valueTitle()
// ---------------------------------------------------------------------------

char const * valueTitle(RenameSchema s)
{
    if (s == RENAME_SEQ)
        return "sequence";
    else
        return "numeric";
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = 0;
    
    // Parse command line.
    FxRenamerOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 1)
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    \t" << options.verbosity << "\n"
                  << "IN           \t" << options.inPath << "\n"
                  << "OUT          \t" << options.outPath << "\n"
                  << "RENAME SCHEMA\t" << valueTitle(options.renameSchema) << "\n";

    // -----------------------------------------------------------------------
    // Open Files.
    // -----------------------------------------------------------------------

    if (options.verbosity >= 1)
        std::cerr << "____OPENING FILES_____________________________________________________________\n"
                  << "\n";

    if (options.verbosity >= 1)
        std::cerr << "INPUT FILE   \t" << options.inPath << "... ";

    seqan::SequenceStream inStream(toCString(options.inPath));
    if (!isGood(inStream))
    {
        std::cerr << "\nERROR: Could not open file " << options.inPath << " for reading!\n";
        return 1;
    }
    if (options.verbosity >= 1)
        std::cerr << "OK\n";

    if (options.verbosity >= 1)
        std::cerr << "OUTPUT FILE  \t" << options.outPath << "... ";
    seqan::SequenceStream outStream(toCString(options.outPath), seqan::SequenceStream::WRITE);
    if (!isGood(outStream))
    {
        std::cerr << "ERROR: Could not open file " << options.outPath << " for writing!\n";
        return 1;
    }
    if (options.verbosity >= 1)
        std::cerr << "OK\n";

    // -----------------------------------------------------------------------
    // Perform Renaming
    // -----------------------------------------------------------------------

    if (options.verbosity >= 1)
        std::cerr << "\n____PERFORMING RENAMING_______________________________________________________\n"
                  << "\n"
                  << "Working...";

    // We also read the sequence into a CharString to this works on IUPAC characters and proteins as well.
    seqan::CharString id, seq, qual;

    for (unsigned i = 1; !atEnd(inStream); ++i)
    {
        if (readRecord(id, seq, qual, inStream) != 0)
        {
            std::cerr << "\nERROR: Error reading record!\n";
            return 1;
        }

        if (options.renameSchema == RENAME_SEQ)
        {
            id = seq;
        }
        else  // RENAME_NUMERIC
        {
            std::stringstream ss;
            ss << i;
            id = ss.str();
        }

        if (writeRecord(outStream, id, seq, qual) != 0)
        {
            std::cerr << "\nERROR: Error writing record!\n";
            return 1;
        }
    }

    std::cerr << "\n";

    if (options.verbosity >= 2)
        std::cerr << "Took " << (sysTime() - startTime) << " s\n";

    return 0;
}
