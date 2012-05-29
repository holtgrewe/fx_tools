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
// Swiss Army Knife tool... "It slices, it dices and it makes the laundry!"
//
// Rewrite of the original sak tool.
// ==========================================================================

#include <sstream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

// --------------------------------------------------------------------------
// Class FxSakOptions
// --------------------------------------------------------------------------

struct FxSakOptions
{
    // Verbosity level.  0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;
    
    // Path to FASTA/FASTQ file.
    seqan::CharString inFastxPath;

    // Path to output file.
    seqan::CharString outPath;

    // Whether or not to print FASTQ to output.
    bool outFastq;

    // Set if one sequence is to be retrieved.
    seqan::String<__uint64> seqIndices;

    // Set if multiple sequences are to be retrieved.
    seqan::String<seqan::Pair<__uint64> > seqIndexRanges;

    // Set if output is to be limited to an infix.
    __uint64 seqInfixBegin;
    __uint64 seqInfixEnd;

    // Whether or not to reverse-complement the result.
    bool reverseComplement;

    // Maximal length of sequence characters to print.
    __uint64 maxLength;

    // Prefix of read names to output if not empty.
    seqan::CharString readPattern;

    FxSakOptions() :
            verbosity(1),
            outFastq(false),
            seqInfixBegin(seqan::maxValue<__uint64>()),
            seqInfixEnd(seqan::maxValue<__uint64>()),
            reverseComplement(false),
            maxLength(seqan::maxValue<__uint64>())
    {}
};

// --------------------------------------------------------------------------
// Function parseRange()
// --------------------------------------------------------------------------

template <typename TNum>
bool parseRange(TNum & beginPos, TNum & endPos, seqan::CharString const & rangeStr)
{
    seqan::Stream<seqan::CharArray<char const *> > stream(begin(rangeStr, seqan::Standard()),
                                                          end(rangeStr, seqan::Standard()));
    seqan::RecordReader<seqan::Stream<seqan::CharArray<char const *> >, seqan::SinglePass<> > reader(stream);

    // Parse out begin position.
    seqan::CharString buffer;
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
    if (!lexicalCast2(beginPos, buffer))
        return false;
    if (beginPos <= 0)
        return false;
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
    if (!lexicalCast2(endPos, buffer))
        return false;
    if (endPos <= 0)
        return false;
    return true;
}

// --------------------------------------------------------------------------
// Function parseArgs()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseArgs(FxSakOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_sak");
    setShortDescription(parser, "Slicing and dicing of FASTA/FASTQ files..");
    setVersion(parser, "0.1");
    setDate(parser, "May 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIIN.fx\\fP");
    addDescription(parser, "\"It slices, it dices and it makes the laundry!\"");

    // The only argument is the input file.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUTFILE, false, "IN"));

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose, log to STDERR."));
    hideOption(parser, "verbose");
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption("o", "out-path", "Path to the resulting file.  If omitted, result is printed to stdout.", seqan::ArgParseArgument::STRING, false, "FASTX"));
    addOption(parser, seqan::ArgParseOption("q", "qual", "Write output as  FASTQ file."));
    addOption(parser, seqan::ArgParseOption("rc", "revcomp", "Reverse-complement output."));
    addOption(parser, seqan::ArgParseOption("l", "max-length", "Maximal number of sequence characters to write out.", seqan::ArgParseArgument::INTEGER, false, "LEN"));

    addSection(parser, "Filter Options");
    addOption(parser, seqan::ArgParseOption("s", "sequence", "Select the given sequence for extraction by 0-based index.", seqan::ArgParseArgument::INTEGER, true, "NUM"));
    addOption(parser, seqan::ArgParseOption("sn", "sequence-name", "Select sequence with name prefix being \\fINAME\\fP.", seqan::ArgParseArgument::STRING, true, "NAME"));
    addOption(parser, seqan::ArgParseOption("ss", "sequences", "Select sequences \\fIfrom\\fP-\\fIto\\fP where \\fIfrom\\fP and \\fIto\\fP are 0-based indices.", seqan::ArgParseArgument::STRING, true, "RANGE"));
    addOption(parser, seqan::ArgParseOption("i", "infix", "Select characters \\fIfrom\\fP-\\fIto\\fP where \\fIfrom\\fP and \\fIto\\fP are 0-based indices.'", seqan::ArgParseArgument::STRING, true, "RANGE"));

    addTextSection(parser, "Usage Examples");
    addListItem(parser, "\\fBfx_sak\\fP \\fB-s\\fP \\fI10\\fP \\fIIN.fa\\fP", "Cut out 11th sequence from \\fIIN.fa\\fP and write to stdout as FASTA.");
    addListItem(parser, "\\fBfx_sak\\fP \\fB-q\\fP \\fB-ss\\fP \\fI10-12\\fP \\fB-ss\\fP \\fI100-200\\fP \\fIIN.fq\\fP", "Cut out 11th up to and including 12th and 101th up to and including 199th sequence from \\fIIN.fq\\fP and write to stdout as FASTQ.");

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        getArgumentValue(options.inFastxPath, parser, 0);

        options.outFastq = isSet(parser, "qual");

        if (isSet(parser, "out-path"))
            getOptionValue(options.outPath, parser, "out-path");

        if (isSet(parser, "verbose"))
            options.verbosity = 2;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 3;

        if (isSet(parser, "sequence"))
        {
            std::vector<std::string> sequenceIds = getOptionValues(parser, "sequence");
            for (unsigned i = 0; i < seqan::length(sequenceIds); ++i)
            {
                unsigned idx = 0;
                if (!seqan::lexicalCast2(idx, sequenceIds[i]))
                {
                    std::cerr << "ERROR: Invalid sequence index " << sequenceIds[i] << "\n";
                    return seqan::ArgumentParser::PARSE_ERROR;
                }
                appendValue(options.seqIndices, idx);
            }
        }

        if (isSet(parser, "sequences"))
        {
            std::vector<std::string> sequenceRanges = getOptionValues(parser, "sequences");
            seqan::CharString buffer;
            for (unsigned i = 0; i < seqan::length(sequenceRanges); ++i)
            {
                seqan::Pair<__uint64> range;
                if (!parseRange(range.i1, range.i2, sequenceRanges[i]))
                {
                    std::cerr << "ERROR: Invalid range " << sequenceRanges[i] << "\n";
                    return seqan::ArgumentParser::PARSE_ERROR;
                }
                appendValue(options.seqIndexRanges, range);
            }
        }

        if (isSet(parser, "infix"))
        {
            seqan::CharString buffer;
            getOptionValue(buffer, parser, "infix");
            if (!parseRange(options.seqInfixBegin, options.seqInfixEnd, buffer))
            {
                std::cerr << "ERROR: Invalid range " << buffer << "\n";
                return seqan::ArgumentParser::PARSE_ERROR;
            }
        }

        options.reverseComplement = isSet(parser, "revcomp");
        
        if (isSet(parser, "max-length"))
            getOptionValue(options.maxLength, parser, "max-length");

        if (isSet(parser, "sequence-name"))
            getOptionValue(options.readPattern, parser, "sequence-name");
    }

    return res;
}

// ---------------------------------------------------------------------------
// Function yesNo()
// ---------------------------------------------------------------------------

char const * yesNo(bool b)
{
    if (b)
        return "YES";
    else
        return "NO";
}

// ---------------------------------------------------------------------------
// Function main()
// ---------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    double startTime = 0;
    
    // Parse command line.
    FxSakOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // -----------------------------------------------------------------------
    // Show options.
    // -----------------------------------------------------------------------
    if (options.verbosity >= 2)
    {
        std::cerr << "____OPTIONS___________________________________________________________________\n"
                  << "\n"
                  << "VERBOSITY    " << options.verbosity << "\n"
                  << "IN           " << options.inFastxPath << "\n"
                  << "OUT          " << options.outPath << "\n"
                  << "FASTQ OUT    " << yesNo(options.outFastq) << "\n"
                  << "INFIX BEGIN  " << options.seqInfixBegin << "\n"
                  << "INFIX END    " << options.seqInfixEnd << "\n"
                  << "MAX LEN      " << options.maxLength << "\n"
                  << "READ PATTERN " << options.readPattern << "\n"
                  << "REVCOMP      " << yesNo(options.reverseComplement) << "\n"
                  << "SEQUENCES\n";
        for (unsigned i = 0; i < length(options.seqIndices); ++i)
            std::cerr << "  SEQ  " << options.seqIndices[i] << "\n";
        for (unsigned i = 0; i < length(options.seqIndexRanges); ++i)
            std::cerr << "  SEQS " << options.seqIndexRanges[i].i1 << "-" << options.seqIndexRanges[i].i2 << "\n";
    }

    // -----------------------------------------------------------------------
    // Open File.
    // -----------------------------------------------------------------------
    std::ostream * outPtr = & std::cout;
    std::fstream inStream;
    if (!empty(options.inFastxPath))
    {
        inStream.open(options.inFastxPath, std::ios::binary | std::ios::in);
        if (!inStream.good())
        {
            std::cerr << "ERROR: Could not open input file " << options.inFastxPath << "\n";
            return 1;
        }
        outPtr = &inStream;
    }

    // Compute index of last sequence to write if any.
    __uint64 endIdx = maxValue<__uint64>();
    SEQAN_FAIL("Write me!");

    // -----------------------------------------------------------------------
    // Read and Write Filtered.
    // -----------------------------------------------------------------------
    startTime = sysTime();
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);
    seqan::AutoSeqStreamFormat tagSelector;
    bool b = checkStreamFormat(reader, tagSelector);
    if (tagSelector.tagId != 1 && tagSelector.tagId != 2)
    {
        std::cerr << "ERROR: Could not determine input format!\n";
        return 1;
    }

    unsigned idx = 0;
    __uint64 charsWritten = 0;
    seqan::CharString id;
    seqan::CharString seq;
    seqan::CharString quals;
    while (!atEnd(reader) && charsWritten < options.maxLength && idx < endIdx)
    {
        if (tagSelector.tagId == 1)
        {
            // FASTA.
            if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
            {
                std::cerr << "ERROR: Reading record!\n";
                return 1;
            }
            if (options.outFastq)
                resize(quals, length(seq), 'I');
        }
        else
        {
            // FASTQ
            if (readRecord(id, seq, quals, reader, seqan::Fastq()) != 0)
            {
                std::cerr << "ERROR: Reading record!\n";
                return 1;
            }
        }

        // Check whether to write out sequence.
        bool writeOut = false;
        // One of options.seqIndices.
        for (unsigned i = 0; i < length(options.seqIndices); ++i)
        {
            if (options.seqIndices[i] == idx)
            {
                writeOut = true;
                break;
            }
        }
        // One of options.seqIndexRanges.
        if (!writeOut)
        {
            for (unsigned i = 0; i < length(options.seqIndexRanges); ++i)
            {
                if (idx >= options.seqIndexRanges[i].i1 && idx < options.seqIndexRanges[i].i2)
                {
                    writeOut = true;
                    break;
                }
            }
        }
        // Name pattern matches.
        if (!writeOut && !empty(options.readPattern))
        {
            unsigned l = length(options.readPattern);
            if (l > length(id))
                l = length(id);
            if (prefix(id, l) == prefix(options.readPattern, l))
                writeOut = true;
        }

        // Write out if we want this.
        if (writeOut)
        {
            if (options.outFastq)
            {
                if (writeRecord(*outPtr, id, seq, quals, Fastq()) != 0)
                {
                    std::cerr << "ERROR: Writing record!\n";
                    return 1;
                }
            }
            else
            {
                if (writeRecord(*outPtr, id, seq, Fasta()) != 0)
                {
                    std::cerr << "ERROR: Writing record!\n";
                    return 1;
                }
            }
        }

        // Advance counter idx.
        idx += 1;
    }
    // b is true if any format was detected successfully.
if (tagSelector.tagId == 1)
    std::cerr << "Detected FASTA." << std::endl;
else if (tagSelector.tagId == 2)
    std::cerr << "Detected FASTQ." << std::endl;
else
    std::cerr << "Unknown file format!" << std::endl;
    if (options.verbosity >= 2)
        std::cerr << "Took " << (sysTime() - startTime) << " s\n";

    return 0;
}
