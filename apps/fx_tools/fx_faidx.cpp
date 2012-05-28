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
// FASTA indexing and indexed access to regions in FASTA files.
//
// This is the equivalent of the "samtools faidx" command.
// ==========================================================================

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

// TODO(holtgrew): This should go from rabema app into core...
#include "../../../../core/apps/rabema/fai_index.h"

// --------------------------------------------------------------------------
// Class FxFaidxOptions
// --------------------------------------------------------------------------

struct FxFaidxOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;
    
    // Path to FASTA file.
    seqan::CharString inFastaPath;

    // Path to .fai file.
    seqan::CharString inFaiPath;

    // Path to results, stdout if empty.
    seqan::CharString outFastaPath;

    // List of regions to retrieve.
    seqan::String<seqan::CharString> regions;

    FxFaidxOptions() : verbosity(1)
    {}
};

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxFaidxOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_faidx");
    setShortDescription(parser, "Indexing FASTA and indexed FASTA access.");
    setVersion(parser, "0.1");
    setDate(parser, "May 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fB-f\\fP \\fIFASTA\\fP] [\\fB-r\\fP \\fIREGION\\fP]+");
    addDescription(parser, "Equivalent program to samtools faidx.");

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "FASTA / FAIDX Files");
    addOption(parser, seqan::ArgParseOption("f", "fasta-file", "Path to the FASTA file.", seqan::ArgParseArgument::STRING, false, "FASTA"));
    setRequired(parser, "fasta-file");
    addOption(parser, seqan::ArgParseOption("i", "index-file", "Path to the .fai index file.  Defaults to FASTA.fai", seqan::ArgParseArgument::STRING, false, "FASTA"));
    addOption(parser, seqan::ArgParseOption("o", "out-file", "Path to the resulting file.  If omitted, result is printed to stdout.", seqan::ArgParseArgument::STRING, false, "FASTA"));

    addSection(parser, "Regions");
    addOption(parser, seqan::ArgParseOption("r", "region", "Region to retrieve from FASTA file.  You can specify multiple regions with multiple \\fB-r\\fP \\fIREGION\\fP.  Note that regions are one-based, see below for detailed information about the format.", seqan::ArgParseArgument::STRING, true, "REGION"));

    addTextSection(parser, "Regions");
    addText(parser,
            "Regions can be specified in the formats \\fICHR\\fP, \\fICHR\\fP:\\fISTART\\fP, \\fICHR\\fP:\\fISTART\\fP:\\fIEND\\fP.  \\fICHR\\fP is the id of the reference sequence in the FASTA file, \\fISTART\\fP and \\fIEND\\fP are the start end end positions of the region.  These positions are one-based.");
    addTextSection(parser, "Region Examples");
    addListItem(parser, "\\fIchr1\\fP", "All of the sequence with the identifier \"chr1\".");
    addListItem(parser, "\\fIchrX\\fP:\\fI1,000\\fP", "The characters in the X chromsome, starting with the 1,000th base.");
    addListItem(parser, "\\fIchr2\\fP:\\fI1,500,000\\fP-\\fI2,000,000\\fP", "The character 1,500,000 up to and including character 2,000,000 in the same chromosome.");

    addTextSection(parser, "Usage Examples");
    addListItem(parser, "\\fBfx_faidx\\fP \\fB-f\\fP \\fIREF.fa\\fP", "Create index for file \\fIREF.fa\\fP, index is written to \\fIREF.fa.fai\\fP");
    addListItem(parser, "\\fBfx_faidx\\fP \\fB-f\\fP \\fIREF.fa\\fP \\fB-i\\fP \\fIINDEX.fai\\fP", "Create index for file \\fIREF.fa\\fP, index is written to \\fIINDEX.fai\\fP");
    addListItem(parser, "\\fBfx_faidx\\fP \\fB-f\\fP \\fIREF.fa\\fP \\fB-r\\fP \\fIchr1\\fP", "Retrieve sequence named \"chr1\" from file \\fIREF.fa\\fP using the index with the default name \\fIREF.fa.fai\\fP.  The index file name is created if it does not exist.");
    addListItem(parser, "\\fBfx_faidx\\fP \\fB-f\\fP \\fIREF.fa\\fP \\fB-r\\fP \\fIchr1:100-1100\\fP", "Retrieve characters 100 to 1,100 from the sequence named \"chr1\" from file \\fIREF.fa\\fP using the index with the default name \\fIREF.fa.fai\\fP.");
    addListItem(parser, "\\fBfx_faidx\\fP \\fB-f\\fP \\fIREF.fa\\fP \\fB-r\\fP \\fIchr1:100-1100\\fP \\fB-r\\fP \\fIchr2:2,000\\fP", "Retrieve characters 100-1,000 from \"chr1\" and all characters from 2,000 of \"chr2\".");
    
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getOptionValue(options.inFastaPath, parser, "fasta-file");

        // Set default FAI file name.
        options.inFaiPath = options.inFastaPath;
        append(options.inFaiPath, ".fai");
        // Get FAI file name from parser if set.
        if (isSet(parser, "index-file"))
            getOptionValue(options.inFaiPath, parser, "index-file");

        if (isSet(parser, "region"))
            options.regions = getOptionValues(parser, "region");

        if (isSet(parser, "out-file"))
            getOptionValue(options.outFastaPath, parser, "out-file");

        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;
    }

    return res;
}

// ---------------------------------------------------------------------------
// Class Region
// ---------------------------------------------------------------------------

struct Region
{
    // Name of sequence.
    seqan::CharString seqName;
    // Index of sequence in FASTA file.  -1 if not set.
    __int32 seqId;
    // 0-based begin position.  -1 if not set.
    __int32 beginPos;
    // 0-based, C-style end position.  -1 if not set.
    __int32 endPos;

    Region() : seqId(-1), beginPos(-1), endPos(-1)
    {}
};

// ---------------------------------------------------------------------------
// Function parseRegion()
// ---------------------------------------------------------------------------

// Parse regionString and write to region.  region.seqId will not be set but
// region.seqName will be.  Return true on success.

bool parseRegion(Region & region, seqan::CharString const & regionString)
{
    seqan::Stream<seqan::CharArray<char const *> > stream(begin(regionString, seqan::Standard()),
                                                          end(regionString, seqan::Standard()));
    seqan::RecordReader<seqan::Stream<seqan::CharArray<char const *> >, seqan::SinglePass<> > reader(stream);

    // Parse out sequence name.
    seqan::CharString buffer;
    int res = readUntilChar(buffer, reader, ':');
    if (res != 0 && res != seqan::EOF_BEFORE_SUCCESS)
        return 1;  // Parse error.
    region.seqName = buffer;
    if (atEnd(reader))
        return true;  // Done after parsing the sequence name.

    goNext(reader);  // Skip ':'.

    // Parse out begin position.
    clear(buffer);
    while (!atEnd(reader) && value(reader) != '-')
    {
        if (!isdigit(value(reader)) && value(reader) != ',')
            return false;  // Error parsing.
        if (isdigit(value(reader)))
            appendValue(buffer, value(reader));
        goNext(reader);
    }
    if (empty(buffer))
        return false;
    if (!lexicalCast2(region.beginPos, buffer))
        return false;
    if (region.beginPos <= 0)
        return false;
    region.beginPos -= 1;  // Adjust to 0-based.
    if (atEnd(reader))
        return true;
    goNext(reader);  // Skip '-'.

    // Parse out end position.
    clear(buffer);
    while (!atEnd(reader))
    {
        if (!isdigit(value(reader)) && value(reader) != ',')
            return false;  // Error parsing.
        if (isdigit(value(reader)))
            appendValue(buffer, value(reader));
        goNext(reader);
    }
    if (empty(buffer))
        return false;
    if (!lexicalCast2(region.endPos, buffer))
        return false;
    if (region.endPos <= 0)
        return false;
    region.endPos -= 1;  // Adjust to 0-based.

    return atEnd(reader);
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = 0;
    
    // Parse command line.
    FxFaidxOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // ---------------------------------------------------------------------------
    // Index I/O
    // ---------------------------------------------------------------------------

    // Load index, create if necessary.
    startTime = sysTime();
    seqan::FaiIndex faiIndex;
    if (load(faiIndex, toCString(options.inFastaPath), toCString(options.inFaiPath)) != 0)
    {
        if (options.verbosity >= 2)
            std::cerr << "Building Index        " << options.inFaiPath << " ...";
        if (buildIndex(toCString(options.inFastaPath), toCString(options.inFaiPath), seqan::Fai()) != 0)
        {
            std::cerr << "Could not build FAI index at " << options.inFaiPath
                      << " for FASTA file " << options.inFastaPath << "\n";
            return 1;
        }
        if (load(faiIndex, toCString(options.inFastaPath), toCString(options.inFaiPath)) != 0)
        {
            std::cerr << "Could not load FAI index we just build.\n";
            return 1;
        }
    }
    if (options.verbosity >= 3)
        std::cerr << "Took " << (startTime - sysTime()) << " s\n";

    // ---------------------------------------------------------------------------
    // Parse and Fetch Regions.
    // ---------------------------------------------------------------------------

    if (empty(options.regions))
        return 0;

    // Parse out regions.
    seqan::String<Region> regions;
    for (unsigned i = 0; i < length(options.regions); ++i)
    {
        Region region;
        if (!parseRegion(region, options.regions[i]))
        {
            std::cerr << "Could not parse region " << options.regions[i] << "\n";
            return 1;
        }
        unsigned seqId;
        if (!getIdByName(faiIndex, region.seqName, seqId))
        {
            std::cerr << "Unknown sequence for region " << options.regions[i] << "\n";
            return 1;
        }
        region.seqId = seqId;
        if (region.seqId < 0 || (unsigned)region.seqId >= length(faiIndex.indexEntryStore))
        {
            std::cerr << "Invalid region " << options.regions[i] << "\n";
            return 1;
        }
        appendValue(regions, region);
    }

    // Open output file.
    std::ostream * outPtr = &std::cout;
    std::ofstream outF;
    if (!empty(options.outFastaPath))
    {
        outF.open(toCString(options.outFastaPath), std::ios::binary | std::ios::out);
        if (!outF.good())
        {
            std::cerr << "Could not open output file " << options.outFastaPath << "\n";
            return 1;
        }
    }

    // Retrieve output infixes and write to result.
    for (unsigned i = 0; i < length(regions); ++i)
    {
        Region const & region = regions[i];
        seqan::CharString id = options.regions[i];
        seqan::Dna5String seq;
        unsigned beginPos = 0;
        if (region.beginPos > 0)
            beginPos = region.beginPos;
        unsigned endPos = sequenceLength(faiIndex, (unsigned)region.seqId);
        if (region.endPos > 0 && (unsigned)region.endPos < endPos)
            endPos = region.endPos;
        if (beginPos > endPos)
            endPos = beginPos;
        getSequenceInfix(seq, faiIndex, region.seqId, beginPos, endPos);
        if (writeRecord(*outPtr, id, seq, seqan::Fasta()) != 0)
        {
            std::cerr << "Could not write infix for region " << options.regions[i] << " to output.\n";
            return 1;
        }
    }

    return 0;
}
