// ==========================================================================
//                                 FX Tools
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

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The input file name is a string.
    seqan::CharString inFilename;

    // The out file name is an out file.
    seqan::CharString outFilename;

    AppOptions() :
        verbosity(1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("fx_read_fastq");
    // Set short description, version, and date.
    setShortDescription(parser, "FASTQ File Provider.");
    setVersion(parser, "0.1");
    setDate(parser, "August 2012");

    // Define usage line and long description.
    addUsageLine(parser, "\\fIINPUT\\fP \\fIOUTPUT\\fP");
    addDescription(parser, "Read a sequence file and provide it as a FASTQ file.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "INPUT"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUTFILE, "OUT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    seqan::getArgumentValue(options.inFilename, parser, 0);
    seqan::getArgumentValue(options.outFilename, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Open in and out file.
    seqan::SequenceStream inStream(toCString(options.inFilename));
    if (!isGood(inStream))
    {
        std::cerr << "ERROR: Could not open file " << options.inFilename << " for reading.\n";
        return 1;
    }
    seqan::SequenceStream outStream(toCString(options.outFilename), seqan::SequenceStream::WRITE,
                                    seqan::SequenceStream::FASTQ);
    if (!isGood(outStream))
    {
        std::cerr << "ERROR: Could not open file " << options.outFilename << " for writing.\n";
        return 1;
    }

    // Copy over the sequence file record wise into a FASTQ file.
    seqan::CharString id;
    seqan::Dna5String seq;
    while (!atEnd(inStream))
    {
        if (readRecord(id, seq, inStream) != 0)
        {
            std::cerr << "ERROR: Error reading from " << options.inFilename << ".\n";
            return 1;
        }
        if (writeRecord(outStream, id, seq) != 0)
        {
            std::cerr << "ERROR: Error writing to " << options.outFilename << ".\n";
            return 1;
        }
    }

    return 0;
}
