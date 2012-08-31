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
// Conversion for FASTA and FASTQ files.  FASTA to FASTQ conversion is less
// interesting than FASTQ than FASTA conversion, of course, but FASTQ to
// FASTQ conversion is useful when converting qualities.
// ==========================================================================

// Reference: Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM. The Sanger FASTQ
// file format for sequences with quality scores, and the Solexa/Illumina FASTQ
// variants. Nucl. Acids Res. (2010) 38(6): 1767-1771.

// TODO(holtgrew): Rename sanger, solexa, illumina to fastq-sanger, fastq-solexa, fastq-illumina?

#include <cmath>
#include <sstream>
#include <cstdio>

#include <zlib.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

// ===========================================================================
// Argument Parsing
// ===========================================================================

struct FxConvertOptions
{
    // Flag whether to rename to numbers.
    bool renameToNumbers;
    // Flag whether to keep the sequence with Ns.
    bool keepNs;
    // Verbosity: 0 is quiet (default), 1 prints report, 2 prints logging.
    int verbosity;
    // Flag whether to compress output with gzip.
    bool gzip;
    // Flag whether to only guess format and quality scale and exit.
    bool guessFormat;
    // Path to input.
    seqan::CharString inPath;
    // Path to output.
    seqan::CharString outPath;
    // Buffer size.  Cannot be set from the outside at the moment.
    unsigned bufferSize;

    enum Format
    {
        AUTO,
        FASTA,
        FASTQ_ILLUMINA,
        FASTQ_SANGER,
        FASTQ_SOLEXA
    };

    Format sourceFormat;
    Format targetFormat;

    FxConvertOptions() : renameToNumbers(false), keepNs(false), verbosity(0), gzip(false), guessFormat(false),
                         bufferSize(4096), sourceFormat(AUTO), targetFormat(FASTA)
    {}
};

// Parse arguments and store them in options.

seqan::ArgumentParser::ParseResult
parseArgs(FxConvertOptions & options,
          int argc,
          char const ** argv)
{
    seqan::ArgumentParser parser("fx_convert");
    setShortDescription(parser, "Sequence File Conversion");
    setVersion(parser, "0.1");
    setDate(parser, "Apr 2012");
    
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fB-i\\fP \\fIINFILE\\fP] [\\fB-o\\fP \\fIOUTFILE\\fP]");
    addDescription(parser, "This is a replacement for fastq_to_fasta from the FASTX toolkit with some extensions.");

    // TODO(holtgrew): I want a custom help text!
    // addOption(parser, seqan::ArgParseOption("h", "help", "This helpful screen."));
    addOption(parser, seqan::ArgParseOption("v", "verbose",
                                            "Verbose - report number of sequences.  If -o is specified, the report "
                                            "will be printed to STDOUT, to STDERR otherwise."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose, log to STDERR."));
    hideOption(parser, "very-verbose");

    addSection(parser, "Filter Related");
    addOption(parser, seqan::ArgParseOption("r", "rename-to-numbers", "Rename sequence identifiers to numbers."));
    addOption(parser, seqan::ArgParseOption("n", "keep-with-ns", "Keep sequences with unknown (N) nucleotides.  "
                                            "Default is to discard such sequences."));

    addSection(parser, "I/O Related");
    addOption(parser, seqan::ArgParseOption("z", "gzip", "Compress output with GZIP."));
    addOption(parser, seqan::ArgParseOption("i", "in-file", "Input file name.", seqan::ArgParseArgument::INPUTFILE));
    setValidValues(parser, "in-file", "fastq fq fasta fa");
    addOption(parser, seqan::ArgParseOption("o", "out-file", "Output file name.", seqan::ArgParseArgument::OUTPUTFILE));
    setValidValues(parser, "out-file", "fastq fq fasta fa");

    addSection(parser, "Quality Related");
    addOption(parser, seqan::ArgParseOption("g", "guess-format", "Guess format and quality scale and exit."));
    addOption(parser, seqan::ArgParseOption("s", "source-format",
                                            "Source quality scale for FASTQ, 'fasta', see Quality Remarks. "
                                            "One of {fasta, sanger, solexa, illumina}.  By default, the input "
                                            "format is detected automatically.", seqan::ArgParseArgument::STRING));
    setValidValues(parser, "source-format", "fasta sanger solexa illumina");
    addOption(parser, seqan::ArgParseOption("t", "target-format", "Target quality scale for FASTQ or 'fasta' (default), "
                                            "see Quality Remarks. One of {fasta, sanger, solexa, illumina}.",
                                            seqan::ArgParseArgument::STRING));
    setValidValues(parser, "target-format", "fasta sanger solexa illumina");

    addTextSection(parser, "Quality Remarks");
    addText(parser,
            "There are four variants for storing qualities in FASTQ files: (1) Sanger-style, storing PHRED "
            "qualities, (2) Solexa-style, (3) Illumina-style.  For more information see the Wikipedia page "
            "on the FASTQ format:");
    addText(parser,
            "If the input is FASTA then the output will be FASTA as well, any values to \\fB-s\\fP and "
            "\\fB-t\\fP are ignored.");
    addText(parser, "http://en.wikipedia.org/wiki/FASTQ_format");
              
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBfx_convert\\fP \\fB-g\\fP < \\fIIN.fq\\fP",
                "Read file \\fIIN.fq\\fP, guess format, print it to stdout and exit.");
    addListItem(parser, "\\fBfx_convert\\fP \\fB-i\\fP \\fIIN.fq\\fP \\fB-o\\fP \\fIOUT.fa\\fP",
                "Read file \\fIIN.fq\\fP, write out as FASTA to \\fIOUT.fa\\fP.");
    addListItem(parser,
                "\\fBfx_convert\\fP \\fB-i\\fP \\fIIN.fq\\fP \\fB-s\\fP \\fIsolexa\\fP \\fB-o\\fP \\fIOUT.fq\\fP "
                " \\fB-s\\fP \\fIsanger\\fP",
                "Read file \\fIIN.fq\\fP with the hint that it has Solexa qualities and write out to \\fIOUT.fq\\fP "
                "using Sanger qualities.");
    
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res == seqan::ArgumentParser::PARSE_OK)
    {
        options.renameToNumbers = isSet(parser, "rename-to-numbers");
        options.keepNs = isSet(parser, "keep-with-ns");
        if (isSet(parser, "verbose"))
            options.verbosity = 1;
        if (isSet(parser, "very-verbose"))
            options.verbosity = 2;
        options.guessFormat = isSet(parser, "guess-format");
        options.gzip = isSet(parser, "gzip");

        if (isSet(parser, "source-format"))
        {
            seqan::CharString tmp;
            getOptionValue(tmp, parser, "source-format");
            if (tmp == "fasta")
                options.sourceFormat = FxConvertOptions::FASTA;
            else if (tmp == "illumina")
                options.sourceFormat = FxConvertOptions::FASTQ_ILLUMINA;
            else if (tmp == "sanger")
                options.sourceFormat = FxConvertOptions::FASTQ_SANGER;
            else if (tmp == "solexa")
                options.sourceFormat = FxConvertOptions::FASTQ_SOLEXA;
            else
                SEQAN_FAIL("Invalid valid for --source-format: %s!", toCString(tmp));
        }
        
        if (isSet(parser, "target-format"))
        {
            seqan::CharString tmp;
            getOptionValue(tmp, parser, "target-format");
            if (tmp == "fasta")
                options.targetFormat = FxConvertOptions::FASTA;
            else if (tmp == "illumina")
                options.targetFormat = FxConvertOptions::FASTQ_ILLUMINA;
            else if (tmp == "sanger")
                options.targetFormat = FxConvertOptions::FASTQ_SANGER;
            else if (tmp == "solexa")
                options.targetFormat = FxConvertOptions::FASTQ_SOLEXA;
            else
                SEQAN_FAIL("Invalid valid for --target-format: %s!", toCString(tmp));
        }
        
        getOptionValue(options.inPath, parser, "in-file");
        getOptionValue(options.outPath, parser, "out-file");
    }

    return res;
}

// ===========================================================================
// Quality Guessing
// ===========================================================================

// The following is reproduced from Wikipedia (http://en.wikipedia.org/wiki/FASTQ_format).  We try to guess the quality
// by looking at the qualities and exclude possibilities.  We fold Illumina 1.8+ and Sanger into the Sanger format and Illumina 1.3+ and Illumina 1.5+ into Illumina.
//
//  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
//  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
//  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
//  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
//  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
//  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
//  |                         |    |        |                              |                     |
// 33                        59   64       73                            104                   126
//
// S - Sanger        Phred+33,  raw reads typically (0, 40)
// X - Solexa        Solexa+64, raw reads typically (-5, 40)
// I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
// J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
//    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
//    (Note: See discussion above).
// L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

struct QualityFormatGuess
{
    // The flags give whether the given format is possible.
    bool sanger;
    bool solexa;
    bool illumina;

    enum BestGuess
    {
        NONE,
        SANGER,
        SOLEXA,
        ILLUMINA
    };

    QualityFormatGuess() : sanger(true), solexa(true), illumina(true)
    {}
};

// Return best guess from quality format guess.

QualityFormatGuess::BestGuess bestGuess(QualityFormatGuess const & guess)
{
    // We can only give a best guess if there is exactly one remaining possibility.
    if (guess.sanger + guess.solexa + guess.illumina != 1)
        return QualityFormatGuess::NONE;
    if (guess.sanger)
        return QualityFormatGuess::SANGER;
    if (guess.solexa)
        return QualityFormatGuess::SOLEXA;
    if (guess.illumina)
        return QualityFormatGuess::ILLUMINA;
    return QualityFormatGuess::NONE;  // Will never reach here.
}

// Return false on errors, true on successful update, regardless of whether any format is still possible.

bool updateQualityFormatGuess(QualityFormatGuess & guess, seqan::CharString const & quals)
{
    typedef seqan::Iterator<seqan::CharString const, seqan::Rooted>::Type TIter;

    for (TIter it = begin(quals, seqan::Rooted()); !atEnd(it); goNext(it))
    {
        char c = *it;
        if (c < 33 && c > 104)
        {
            std::cerr << "ERROR: Invalid quality " << (int)c << "!\n";
            return false;
        }
        if (c > 74)  // Sanger allows <= 73, but Illumina 1.8 allows <= 74.
            guess.sanger = false;
        if (c < 59)
            guess.solexa = false;
        if (c < 64)
            guess.illumina = false;
    }
    return true;
}

// ===========================================================================
// Quality Conversion
// ===========================================================================

// Given an input and an output format, filles a quality conversion array.  The array has 256 entries and maps source
// quality chars to target quality chars.

void buildConversionTable(char * table,
                          QualityFormatGuess::BestGuess source,
                          QualityFormatGuess::BestGuess target)
{
    // Get minimum and maximum source value.
    unsigned from = 0, to = 0;
    if (source == QualityFormatGuess::SANGER)
    {
        from = 33;
        to = 126;
    }
    else if (source == QualityFormatGuess::SOLEXA)
    {
        from = 59;
        to = 126;
    }
    else if (source == QualityFormatGuess::ILLUMINA)
    {
        from = 64;
        to = 126;
    }
    else
    {
        SEQAN_FAIL("ERROR: Cannot build conversion table for unknown source.");
    }

    // Initialize, everything with minimum or maximum value.
    for (unsigned i = 0; i < 256; ++i)
    {
        if (target == QualityFormatGuess::SANGER)
            table[i] = (i <= from) ? 33 : 73;
        else if (target == QualityFormatGuess::SOLEXA)
            table[i] = (i <= from) ? 59 : 104;
        else if (target == QualityFormatGuess::ILLUMINA)
            table[i] = (i <= from) ? 64 : 104;
    }

    for (unsigned i = from; i <= to; ++i)
    {
        // Convert from quality i to error probability.
        double p = 0;
        if (source == QualityFormatGuess::SANGER)
            p = pow(10.0, ((int)i - 33) / -10.0);
        else if (source == QualityFormatGuess::SOLEXA)
            p = 1.0 / (pow(10, ((int)i - 64) / 10.0) + 1);
        else // if (source == QualityFormatGuess::ILLUMINA)
            p = pow(10.0, ((int)i - 64) / -10.0);
        SEQAN_ASSERT_GEQ(p, 0.0);
        SEQAN_ASSERT_LEQ(p, 1.0);

        // Now convert from error probability to quality.
        int q = 0;
        if (target == QualityFormatGuess::SANGER)
        {
            q = -10.0 * log10(p);
            if (q < 0)
                q = 0;
            if (q > 93)
                q = 93;
            q += 33;
        }
        else if (target == QualityFormatGuess::SOLEXA)
        {
            q = -10 * log10(p / (1.0 - p));
            if (q < -5)
                q = -5;
            if (q > 62)
                q = 62;
            q += 64;
        }
        else // if (target == QualityFormatGuess::ILLUMINA)
        {
            q = -10.0 * log10(p);
            if (q < 0)
                q = 62;
            if (q > 0)
                q = 62;
            q += 64;
        }
        SEQAN_ASSERT_GT(q, 0);
        SEQAN_ASSERT_LT(q, 255);
        table[i] = static_cast<char>(q);

        // std::cerr << "i = " << i << "\tq = " << q << "\tp = " << p << '\n';
    }
}

// ===========================================================================
// Main Program
// ===========================================================================

template <typename TOutStream>
int runConvert(TOutStream & out,
               std::ostream & err,
               std::istream & in,
               FxConvertOptions const & options)
{
    typedef seqan::RecordReader<std::istream, seqan::SinglePass<> > TRecordReader;
    TRecordReader reader(in);

    // TODO(holtgrew): Check whether it would be faster to multiplex to formats before calling readRecord.
    
    // Guess format.
    seqan::AutoSeqStreamFormat tagSelector;
    if (!checkStreamFormat(reader, tagSelector))
    {
        err << "ERROR: Cannot determine file format.\n";
        return 1;
    }
    if (options.verbosity >= 2)
    {
        if (tagSelector.tagId == 1)
            err << "File format is FASTA.\n";
        else
            err << "File format is FASTQ.\n";
    }

    // If the format is FASTQ, then guess the quality format.
    if (tagSelector.tagId == 2)
    {
        QualityFormatGuess qualityFormatGuess;
        seqan::CharString id, qual;
        seqan::Dna5String seq;
        QualityFormatGuess::BestGuess formatGuess;
            
        // Guess input quality format if not specified.
        if (options.sourceFormat == FxConvertOptions::AUTO)
        {
            // Create scope for record reader limiting.
            {
                seqan::LimitRecordReaderInScope<std::istream, seqan::SinglePass<> > limiter(reader);
                while (!atEnd(reader))
                {
                    int res = readRecord(id, seq, qual, reader, seqan::Fastq());
                    if (res != 0 && res != seqan::EOF_BEFORE_SUCCESS)
                    {
                        std::cerr << "Error reading input!\n";
                        return 1;
                    }
                    if (!updateQualityFormatGuess(qualityFormatGuess, qual))
                        return 1;
                }
            }

            formatGuess = bestGuess(qualityFormatGuess);
            if (formatGuess == QualityFormatGuess::NONE)
            {
                std::cerr << "ERROR: Could not guess FASTQ quality scale unambiguously!\n";
                if (qualityFormatGuess.sanger)
                    std::cerr << "Could be Sanger.\n";
                if (qualityFormatGuess.solexa)
                    std::cerr << "Could be Solexa.\n";
                if (qualityFormatGuess.illumina)
                    std::cerr << "Could be Illumina.\n";
                return 1;
            }
        }
        else
        {
            switch (options.sourceFormat)
            {
                case FxConvertOptions::FASTQ_ILLUMINA:
                    formatGuess = QualityFormatGuess::ILLUMINA;
                    break;
                case FxConvertOptions::FASTQ_SANGER:
                    formatGuess = QualityFormatGuess::SANGER;
                    break;
                case FxConvertOptions::FASTQ_SOLEXA:
                    formatGuess = QualityFormatGuess::SOLEXA;
                    break;
                default:
                    SEQAN_FAIL("Should never reach here!\n");
            }
        }

        seqan::CharString format;
        switch (formatGuess)
        {
            case QualityFormatGuess::SOLEXA:
                format = "text/x-fastq-solexa";
                break;
            case QualityFormatGuess::ILLUMINA:
                format = "text/x-fastq-illumina";
                break;
            default:
                format = "text/x-fastq-sanger";
        }
        
        // If we only wanted to guess the file format then we are done here.  Otherwise, we only print the file format
        // when the verbosity is high and carry on.
        if (options.guessFormat)
        {
            seqan::streamPut(out, "content-type: ");
            seqan::streamPut(out, format);
            seqan::streamPut(out, '\n');
            return 0;
        }
        if (options.verbosity >= 2)
            err << "Guessed input quality scale to be " << format << "\n";

        // Compute quality conversion table.
        char qualityConversionTable[256];
        QualityFormatGuess::BestGuess outFormat = QualityFormatGuess::NONE;  // FASTA
        switch(options.targetFormat)
        {
            case FxConvertOptions::FASTQ_ILLUMINA:
                outFormat = QualityFormatGuess::ILLUMINA;
                break;
            case FxConvertOptions::FASTQ_SANGER:
                outFormat = QualityFormatGuess::SANGER;
                break;
            case FxConvertOptions::FASTQ_SOLEXA:
                outFormat = QualityFormatGuess::SOLEXA;
                break;
            default:  // FASTA
                outFormat = QualityFormatGuess::NONE;
                break;
        }
        if (outFormat != QualityFormatGuess::NONE)
            buildConversionTable(qualityConversionTable, formatGuess, outFormat);

        // Now, read the whole file and write out.
        std::stringstream ss;
        unsigned num = 1;
        while (!atEnd(reader))
        {
            if (readRecord(id, seq, qual, reader, seqan::Fastq()) != 0)
            {
                std::cerr << "ERROR: Problem reading FASTQ file!\n";
                return 1;
            }

            // Replace sequence identifier if requested.
            if (options.renameToNumbers)
            {
                ss.str("");
                ss.clear();
                ss << num++;
                id = ss.str();
            }

            // TODO(holtgrew): Interpret target format.
            
            if (outFormat == QualityFormatGuess::NONE)
            {
                if (writeRecord(out, id, seq, seqan::Fasta()) != 0)
                {
                    std::cerr << "ERROR: Problem writing FASTA file!\n";
                    return 1;
                }
            }
            else
            {
                // Perform quality scale conversion.
                if (formatGuess != outFormat)
                {
                    // std::stringstream ss;
                    typedef seqan::Iterator<seqan::CharString, seqan::Rooted>::Type TIter;
                    for (TIter it = begin(qual, seqan::Rooted()); !atEnd(it); goNext(it))
                    // {
                        *it = qualityConversionTable[static_cast<__uint8>(*it)];
                    //     ss << "[" << (int)*it << " " << (int)qualityConversionTable[static_cast<__uint8>(*it)] << " " << (char)(int)qualityConversionTable[static_cast<__uint8>(*it)] << "] ";;
                    // }
                    // qual = ss.str();
                }
                
                if (writeRecord(out, id, seq, qual, seqan::Fastq()) != 0)
                {
                    std::cerr << "ERROR: Problem writing FASTQ file!\n";
                    return 1;
                }
            }
        }
    }
    else
    {
        // In the case of FASTA, simply read record by record and write out as FASTA.
        std::stringstream ss;
        unsigned num = 1;
        seqan::CharString id;
        seqan::Dna5String seq;
        while (!atEnd(reader))
        {
            if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
            {
                std::cerr << "ERROR: Problem reading FASTA file!\n";
                return 1;
            }

            // Replace sequence identifier if requested.
            if (options.renameToNumbers)
            {
                ss.str("");
                ss.clear();
                ss << num++;
                id = ss.str();
            }

            if (writeRecord(out, id, seq, seqan::Fasta()) != 0)
            {
                std::cerr << "ERROR: Problem writing FASTA file!\n";
                return 1;
            }
        }
    }

    return 0;
}

// ===========================================================================
// Entry Point
// ===========================================================================

int main(int argc, char const ** argv)
{
    // Parse command line.
    FxConvertOptions options;
    seqan::ArgumentParser::ParseResult res = parseArgs(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;  // 1 on errors, 0 otherwise

    // TODO(holtgrew): Using MMAP strings in case of filename would be faster.
    // TODO(holtgrew): This multiplexing is quite ugly... We could use gzStream for both input and output and attach to 
    // Open input and output stream, forward 
    if (empty(options.inPath))
    {
        if (empty(options.outPath))
        {
            return runConvert(std::cout, std::cerr, std::cin, options);
        }
        else
        {
            if (options.gzip)
            {
                gzFile gzOut = gzopen(seqan::toCString(options.outPath), "wb");
                if (gzOut == NULL)
                {
                    std::cerr << "Could not open " << options.outPath << '\n';
                    return 1;
                }
                seqan::Stream<seqan::GZFile> gzStream(gzOut);
                int res = runConvert(gzStream, std::cerr, std::cin, options);
                gzclose(gzOut);
                return res;
            }
            else
            {
                std::ofstream out(seqan::toCString(options.outPath), std::ios::binary | std::ios::in);
                if (!out.good())
                {
                    std::cerr << "ERROR: Could not open " << options.outPath << '\n';
                    return 1;
                }
                return runConvert(out, std::cerr, std::cin, options);
            }
        }
    }
    else
    {
        std::ifstream in(seqan::toCString(options.inPath), std::ios::binary | std::ios::in);

        if (empty(options.outPath))
        {
            return runConvert(std::cout, std::cerr, in, options);
        }
        else
        {
            if (options.gzip)
            {
                gzFile gzOut = gzopen(seqan::toCString(options.outPath), "wb");
                if (gzOut == NULL)
                {
                    std::cerr << "Could not open " << options.outPath << '\n';
                    return 1;
                }
                seqan::Stream<seqan::GZFile> gzStream(gzOut);
                int res = runConvert(gzStream, std::cerr, in, options);
                gzclose(gzOut);
                return res;
            }
            else
            {
                std::ofstream out(seqan::toCString(options.outPath), std::ios::binary | std::ios::in);
                if (!out.good())
                {
                    std::cerr << "ERROR: Could not open " << options.outPath << '\n';
                    return 1;
                }
                return runConvert(out, std::cerr, in, options);
            }
        }
    }
    
    return 0;
}
