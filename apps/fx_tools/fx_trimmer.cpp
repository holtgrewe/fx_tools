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
// Trimming of read, bases from 5' and 3' end are given as numbers or
// percentages on the command line.
// ==========================================================================

#include <cmath>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// Enum BaseOffsetType
// --------------------------------------------------------------------------

// Whether base offsets are given in nucleotide counts or percentages of read length.

enum BaseOffsetType
{
    OFFSET_COUNT,   // Offsets are counts.
    OFFSET_PERCENT  // Offsets are percentages.
};

// --------------------------------------------------------------------------
// Class FxRenamerOptions
// --------------------------------------------------------------------------

struct FxTrimmerOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;
    
    // Path to FASTA/FASTQ file.
    seqan::CharString inPath;

    // Path to output file.
    seqan::CharString outPath;

    // Whether to use counts or percentages as offsets.
    BaseOffsetType offsetType;

    // Offset from left and right.
    unsigned offsetLeft;
    unsigned offsetRight;

    FxTrimmerOptions() :
            verbosity(1), offsetType(OFFSET_COUNT), offsetLeft(0), offsetRight(0)
    {}
};

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxTrimmerOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_trimmer");
    setShortDescription(parser, "Trim sequences in FASTA and FASTQ files.");
    setVersion(parser, "0.1");
    setDate(parser, "August 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fIOPTIONS\\fP] -i \\fIOUT.fasta\\fP \\fB-o\\fP \\fIOUT.fasta\\fP");
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fIOPTIONS\\fP] -i \\fIOUT.fastq\\fP \\fB-o\\fP \\fIOUT.fastq\\fP");
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

    addSection(parser, "Trimming Options");
    addOption(parser, seqan::ArgParseOption("t", "offset-type", "Select the offset type.", seqan::ArgParseArgument::STRING, "TYPE"));
    setValidValues(parser, "offset-type", "count percentage");
    setDefaultValue(parser, "offset-type", "count");
    addOption(parser, seqan::ArgParseOption("l", "offset-left", "Offset from left end (5'-end).", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "offset-left", "0");
    setDefaultValue(parser, "offset-type", "count");
    addOption(parser, seqan::ArgParseOption("r", "offset-right", "Select the right end (3'-end).", seqan::ArgParseArgument::INTEGER, "INT"));
    setMinValue(parser, "offset-right", "0");
    
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getOptionValue(options.inPath, parser, "in-path");
        getOptionValue(options.outPath, parser, "out-path");

        seqan::CharString offsetType;
        getOptionValue(offsetType, parser, "offset-type");
        if (offsetType == "count")
            options.offsetType = OFFSET_COUNT;
        else
            options.offsetType = OFFSET_PERCENT;

        getOptionValue(options.offsetLeft, parser, "offset-left");
        getOptionValue(options.offsetRight, parser, "offset-right");

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

char const * valueTitle(BaseOffsetType b)
{
    if (b == OFFSET_COUNT)
        return "count";
    else
        return "percent";
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = 0;
    
    // Parse command line.
    FxTrimmerOptions options;
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
                  << "OFFSET TYPE  \t" << valueTitle(options.offsetType) << "\n"
                  << "LEFT OFFSET  \t" << options.offsetLeft << "\n"
                  << "RIGHT OFFSET \t" << options.offsetRight << "\n";

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
        std::cerr << "\n____PERFORMING TRIMMING_______________________________________________________\n"
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

        // Compute offsets.
        unsigned offsetLeft = 0;
        unsigned offsetRight = 0;
        if (options.offsetType == OFFSET_COUNT)
        {
            offsetLeft = options.offsetLeft;
            offsetRight = options.offsetRight;
        }
        else  // RENAME_NUMERIC
        {
            offsetLeft = round(0.01 * options.offsetLeft * length(seq));
            offsetRight = round(0.01 * options.offsetRight * length(seq));
        }

        // Compute start and end position of infix to cut.
        unsigned beginPos = offsetLeft;
        if (offsetLeft > length(seq))
            beginPos = length(seq);
        unsigned endPos = length(seq) - offsetRight;
        if (offsetRight > length(seq))
            endPos = 0;
        if (endPos < beginPos)
            beginPos = endPos;

        // Cut out infix.
        seq = infix(seq, beginPos, endPos);

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
