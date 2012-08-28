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
// Computes coverage and C+G content of a genome given a SAM file,
// window-based.
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

// --------------------------------------------------------------------------
// Class FxSamCoverageOptions
// --------------------------------------------------------------------------

struct FxSamCoverageOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;
    
    // Path to Genome file.
    seqan::CharString inGenomePath;

    // Path to SAM file.
    seqan::CharString inSamPath;

    // Path to output file.
    seqan::CharString outPath;

    // Window size to use for computation.
    __int32 windowSize;

    FxSamCoverageOptions() : verbosity(1), windowSize(10*1000)
    {}
};


// --------------------------------------------------------------------------
// Class BinData
// --------------------------------------------------------------------------

struct BinData
{
    // Number of reads whose alignment starts here.
    unsigned coverage;
    // Length of underlying sequence.
    unsigned length;
    // Percentag C+G.
    double cgContent;
    
    BinData() : coverage(0), length(0), cgContent(0)
    {}
};

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxSamCoverageOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_sam_coverage");
    setShortDescription(parser, "Read Coverage Computation.");
    setVersion(parser, "0.1");
    setDate(parser, "August 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-o\\fP \\fIOUT.tsv\\fP \\fIGENOME.fa\\fP \\fIMAPPING.sam\\fP");
    addDescription(parser, "Compute read coverage and C+G content for a genome.");

    // Two input files: Genome, and mapping.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE));
    setValidValues(parser, 0, "fasta fa");
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE));
    setValidValues(parser, 1, "sam");

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Main Options");
    addOption(parser, seqan::ArgParseOption("w", "window-size", "Set the size of the non-overlapping windows in base pairs.", seqan::ArgParseArgument::INTEGER, "NUM"));
    setDefaultValue(parser, "window-size", "10000");

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-path", "Path to the resulting file.  If omitted, result is printed to stdout.", seqan::ArgParseArgument::OUTPUTFILE, "TSV"));
    setRequired(parser, "out-path");
    setValidValues(parser, "out-path", "sam.coverage.tsv");

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getArgumentValue(options.inGenomePath, parser, 0);
        getArgumentValue(options.inSamPath, parser, 1);
        getOptionValue(options.outPath, parser, "out-path");

        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;
    }

    return res;
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = 0;
    
    // -----------------------------------------------------------------------
    // Parse command line.
    // -----------------------------------------------------------------------
    FxSamCoverageOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 1)
    {
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    " << options.verbosity << "\n"
                  << "GENOME       " << options.inGenomePath << "\n"
                  << "SAM          " << options.inSamPath << "\n"
                  << "OUT          " << options.outPath << "\n"
                  << "WINDOW SIZE  " << options.windowSize << "\n";
    }

    // -----------------------------------------------------------------------
    // Load Genome FAI Index
    // -----------------------------------------------------------------------

    std::cerr << "\n"
              << "___PREPRATION_____________________________________________________________________\n"
              << "\n"
              << "Indexing GENOME file  " << options.inGenomePath << " ...";
    seqan::FaiIndex faiIndex;
    if (build(faiIndex, toCString(options.inGenomePath)) != 0)
    {
        std::cerr << "Could not build FAI index.\n";
        return 1;
    }
    std::cerr << " OK\n";

    // Prepare bins.
    seqan::String<seqan::String<BinData> > bins;
    resize(bins, numSeqs(faiIndex));

    // -----------------------------------------------------------------------
    // Compute C+G content 
    // -----------------------------------------------------------------------

    std::cerr << "\n"
              << "___C+G CONTENT COMPUTATION________________________________________________________\n"
              << "\n";

    for (unsigned i = 0; i < numSeqs(faiIndex); ++i)
    {
        std::cerr << "[" << sequenceName(faiIndex, i) << "] ...";
        unsigned numBins = (sequenceLength(faiIndex, i) + options.windowSize - 1) / options.windowSize;
        resize(bins[i], numBins);
        seqan::Dna5String contigSeq;
        if (readSequence(contigSeq, faiIndex, i) != 0)
        {
            std::cerr << "\nERROR: Could not read sequence " << sequenceName(faiIndex, i) << " from file!\n";
            return 1;
        }

        for (unsigned bin = 0; bin < numBins; ++bin)
        {
            unsigned cgCounter = 0;
            unsigned binSize = 0;
            bins[i][bin].length = options.windowSize;
            if ((bin + 1) * options.windowSize > length(contigSeq))
                bins[i][bin].length = length(contigSeq) - bin * options.windowSize;
            for (unsigned pos = bin * options.windowSize; pos < length(contigSeq) && pos < (bin + 1) * options.windowSize; ++pos, ++binSize)
                cgCounter += (contigSeq[pos] == 'C' || contigSeq[pos] == 'G');
            bins[i][bin].cgContent = 1.0 * cgCounter / binSize;
        }
        std::cerr << "DONE\n";
    }

    // -----------------------------------------------------------------------
    // Compute Coverage
    // -----------------------------------------------------------------------

    std::cerr << "\n"
              << "___COVERAGE COMPUATATION________________________________________________________\n"
              << "\n"
              << "Computing Coverage...";

    seqan::BamStream bamStream(toCString(options.inSamPath));
    if (!isGood(bamStream))
    {
        std::cerr << "Could not open " << options.inSamPath << "!\n";
        return 1;
    }

    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStream))
    {
        if (readRecord(record, bamStream) != 0)
        {
            std::cerr << "ERROR: Could not read record from BAM file!\n";
            return 1;
        }

        if (hasFlagUnmapped(record) || hasFlagSecondary(record) || record.rId == seqan::BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip these records.

        int contigId = 0;
        seqan::CharString const & contigName = nameStore(bamStream.bamIOContext)[record.rId];
        if (!getIdByName(faiIndex, contigName, contigId))
        {
            std::cerr << "ERROR: Alignment to unknown contig " << contigId << "!\n";
            return 1;
        }
        unsigned binNo = record.pos / options.windowSize;
        bins[contigId][binNo].coverage += 1;
    }

    std::cerr << "DONE\n";

    // -----------------------------------------------------------------------
    // Write Output
    // -----------------------------------------------------------------------

    std::ostream * out = &std::cout;
    std::ofstream outFile;
    if (options.outPath != "-")
    {
        outFile.open(toCString(options.outPath), std::ios::binary | std::ios::out);
        if (!outFile.good())
        {
            std::cerr << "ERROR: Could not open output file " << options.outPath << "!\n";
            return 1;
        }
        out = &outFile;
    }

    (*out) << "#BIN\tREF_NAME\tREF_BIN\tBIN_BEGIN\tBIN_LENGTH\tCOVERAGE\tCG_CONTENT\n";
    for (unsigned i = 0, globalBin = 0; i < length(bins); ++i)
    {
        for (unsigned refBin = 0; refBin < length(bins[i]); ++refBin, ++globalBin)
        {
            (*out) << globalBin << '\t'
                   << sequenceName(faiIndex, i) << '\t'
                   << refBin << '\t'
                   << refBin * options.windowSize << '\t'
                   << bins[i][refBin].length << '\t'
                   << bins[i][refBin].coverage << '\t'
                   << bins[i][refBin].cgContent << '\n';
        }
    }

    if (options.verbosity >= 2)
        std::cerr << "Took " << (sysTime() - startTime) << " s\n";

    return 0;
}
