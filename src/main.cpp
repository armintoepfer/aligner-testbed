#include <kallocw.h>
#include <miniwfa.h>
#include <pbcopper/cli2/CLI.h>
#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/utility/Stopwatch.h>
#include <WFAligner.hpp>
#include <boost/algorithm/string.hpp>
#include <pancake/AlignerKSW2.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "pbcopper/data/Cigar.h"
#include "wavefront/wavefront_attributes.h"

using namespace PacBio;

namespace OptionNames {
// clang-format off
const CLI_v2::Option WFA2C {
R"({
    "names" : ["wfa2-c"],
    "description" : "Use WFA2-lib C API",
    "type" : "bool",
    "default" : true
})"
};
const CLI_v2::Option WFA2Cpp {
R"({
    "names" : ["wfa2-cpp"],
    "description" : "Use WFA2-lib C++ API",
    "type" : "bool",
    "default" : false
})"
};
const CLI_v2::Option KSW2 {
R"({
    "names" : ["ksw2"],
    "description" : "Use KSW2",
    "type" : "bool",
    "default" : true
})"
};
const CLI_v2::Option MiniWFA {
R"({
    "names" : ["miniwfa"],
    "description" : "Use miniwfa",
    "type" : "bool",
    "default" : true
})"
};
const CLI_v2::Option Adaptive {
R"({
    "names" : ["wfa-adaptive"],
    "description" : "Use adaptive heuristic for WFA",
    "type" : "bool",
    "default" : false
})"
};
const CLI_v2::Option Rounds {
R"({
    "names" : ["r", "rounds"],
    "description" : "Repeat N rounds and average timings",
    "type" : "int",
    "default" : 1
})"
};
// clang-format on
}  // namespace OptionNames

CLI_v2::Interface CreateCLI()
{
    static const std::string description{"CLR align challenge."};
    const auto version = "0.0.1";
    CLI_v2::Interface i{"cas", description, version};

    Logging::LogConfig logConfig;
    logConfig.Header = "| ";
    logConfig.Delimiter = " | ";
    logConfig.Fields = Logging::LogField::TIMESTAMP | Logging::LogField::LOG_LEVEL;
    logConfig.Level = Logging::LogLevel::INFO;
    i.LogConfig(logConfig);

    const CLI_v2::PositionalArgument InputFile{
        R"({
        "name" : "Pairwise sequences",
        "description" : "TXT.",
        "type" : "file",
        "required" : true
    })"};
    i.AddPositionalArgument(InputFile);
    i.AddOption(OptionNames::Rounds);
    i.AddOption(OptionNames::Adaptive);
    i.AddOptionGroup("Algorithms", {OptionNames::KSW2, OptionNames::WFA2C, OptionNames::WFA2Cpp,
                                    OptionNames::MiniWFA});

    return i;
}

std::vector<std::pair<std::string, std::string>> LoadTargetQueryFromFile(const std::string& str)
{
    std::ifstream is{str};
    if (!is) {
        throw std::runtime_error{"cannot open '" + str + "' for reading"};
    }

    std::vector<std::pair<std::string, std::string>> result;
    for (std::string line; std::getline(is, line);) {
        std::vector<std::string> splitStrs;
        boost::split(splitStrs, line, boost::is_any_of(" "));

        if (splitStrs.size() != 2) {
            throw std::runtime_error{'\'' + line + "' is not a valid target<space>query pair"};
        }

        result.emplace_back(std::move(splitStrs[0]), std::move(splitStrs[1]));
    }

    return result;
}

std::string SplitUnpackedCigar(const std::string& cigarStr)
{
    std::string cigar;
    int num = 0;
    char type = cigarStr[0];
    for (char c : cigarStr) {
        if (isdigit(c)) {
            std::cerr << "Is this a valid unpacked CIGAR? <" << cigarStr << ">?\n";
            exit(1);
        }
        if (c != type) {
            if (type == 'M') {
                type = '=';
            }
            cigar += std::to_string(num) + type;
            type = c;
            num = 0;
        }
        num += 1;
    }
    cigar += std::to_string(num) + type;
    return cigar;
}

int64_t RunMiniWFA(const std::vector<std::pair<std::string, std::string>>& sequences,
                   Pancake::AlignmentParameters& alnP, const int32_t rounds,
                   const Logging::LogLevel logLevel)
{
    int32_t overallCigarLength{0};

    mwf_opt_t opt;
    mwf_opt_init(&opt);
    opt.x = alnP.mismatchPenalty;
    opt.o1 = alnP.gapOpen1;
    opt.e1 = alnP.gapExtend1;
    opt.o2 = alnP.gapOpen2;
    opt.e2 = alnP.gapExtend2;
    opt.flag |= MWF_F_CIGAR;
    void* km = 0;

    km = km_init_w();  // Avoid using ksw2's kalloc..
    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            mwf_rst_t rst;
            mwf_wfa(km, &opt, qry.size(), qry.c_str(), target.size(), target.c_str(), &rst);
            // mwf_assert_cigar(&opt, rst.n_cigar, rst.cigar, qry.size(), target.size(), rst.s);
            std::string cigar;
            for (int32_t j = 0; j < rst.n_cigar; ++j) {
                cigar += std::to_string(rst.cigar[j] >> 4) + "MIDNSHP=XBid"[rst.cigar[j] & 0xf];
            }
            if (logLevel == Logging::LogLevel::DEBUG) {
                PBLOG_DEBUG << cigar;
            }
            overallCigarLength += cigar.size();
            kfree_w(km, rst.cigar);
        }
    }
    km_destroy_w(km);

    return overallCigarLength;
}

int64_t RunWFA2Cpp(const std::vector<std::pair<std::string, std::string>>& sequences,
                   Pancake::AlignmentParameters& alnP, const int32_t rounds,
                   const Logging::LogLevel logLevel)
{
    Utility::Stopwatch timer;
    int32_t overallCigarLength{0};
    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            wfa::WFAlignerGapAffine2Pieces aligner(
                alnP.mismatchPenalty, alnP.gapOpen1, alnP.gapExtend1, alnP.gapOpen2,
                alnP.gapExtend2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
            aligner.alignEnd2End(qry.c_str(), qry.size(), target.c_str(), target.size());
            const std::string cigarCompressed = SplitUnpackedCigar(aligner.getAlignmentCigar());
            overallCigarLength += cigarCompressed.size();
            if (logLevel == Logging::LogLevel::DEBUG) {
                PBLOG_DEBUG << cigarCompressed;
            }
        }
    }
    return overallCigarLength;
}

int64_t RunWFA2C(const std::vector<std::pair<std::string, std::string>>& sequences,
                 Pancake::AlignmentParameters& alnP, const int32_t rounds,
                 const Logging::LogLevel logLevel, const bool adaptive)
{
    int32_t overallCigarLength{0};

    auto attributes = wavefront_aligner_attr_default;
    attributes.memory_mode = wavefront_memory_high;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = alnP.mismatchPenalty;
    attributes.affine2p_penalties.gap_opening1 = alnP.gapOpen1;
    attributes.affine2p_penalties.gap_extension1 = alnP.gapExtend1;
    attributes.affine2p_penalties.gap_opening2 = alnP.gapOpen2;
    attributes.affine2p_penalties.gap_extension2 = alnP.gapExtend2;
    attributes.alignment_scope = compute_alignment;
    auto wf_aligner = wavefront_aligner_new(&attributes);
    if (adaptive) {
        wavefront_aligner_set_heuristic_wfadaptive(wf_aligner, 10, 50, 1);
    } else {
        wavefront_aligner_set_heuristic_none(wf_aligner);
    }
    wavefront_aligner_set_alignment_end_to_end(wf_aligner);

    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            wavefront_align(wf_aligner, qry.c_str(), qry.size(), target.c_str(), target.size());
            char* buffer = wf_aligner->cigar.operations + wf_aligner->cigar.begin_offset;
            int buf_len = wf_aligner->cigar.end_offset - wf_aligner->cigar.begin_offset;
            const std::string cigarLinear(buffer, buf_len);
            const std::string cigarCompressed = SplitUnpackedCigar(cigarLinear);
            overallCigarLength += cigarCompressed.size();
            if (logLevel == Logging::LogLevel::DEBUG) {
                PBLOG_DEBUG << cigarCompressed;
            }
        }
    }
    wavefront_aligner_delete(wf_aligner);
    return overallCigarLength;
}

int64_t RunKSW2(const std::vector<std::pair<std::string, std::string>>& sequences,
                Pancake::AlignmentParameters& alnP, const int32_t rounds,
                const Logging::LogLevel logLevel)
{
    int32_t overallCigarLength{0};

    auto ksw2Aligner = Pancake::CreateAlignerKSW2(alnP);
    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            const auto res = ksw2Aligner->Global(target, qry);
            overallCigarLength += res.cigar.ToStdString().size();
            if (logLevel == Logging::LogLevel::DEBUG) {
                PBLOG_DEBUG << res.cigar.ToStdString();
            }
        }
    }

    return overallCigarLength;
}

int RunnerSubroutine(const CLI_v2::Results& options)
{
    const std::vector<std::string> files = options.PositionalArguments();
    const auto sequences = LoadTargetQueryFromFile(files.front());

    const int32_t seqPairs = sequences.size();
    PBLOG_INFO << "Number of sequence pairs : " << seqPairs;

    Pancake::AlignmentParameters alnP;
    const int32_t rounds = options[OptionNames::Rounds];
    const auto loglevel{options.LogLevel()};

    int64_t cl{0};
    if (options[OptionNames::MiniWFA]) {
        Utility::Stopwatch timer;
        cl += RunMiniWFA(sequences, alnP, rounds, loglevel);
        PBLOG_INFO << "miniwfa time  : "
                   << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() /
                                                                 rounds / sequences.size());
    }
    if (options[OptionNames::WFA2C]) {
        Utility::Stopwatch timer;
        cl += RunWFA2C(sequences, alnP, rounds, loglevel, options[OptionNames::Adaptive]);
        PBLOG_INFO << "WFA2 C time   : "
                   << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() /
                                                                 rounds / sequences.size());
    }
    if (options[OptionNames::WFA2Cpp]) {
        Utility::Stopwatch timer;
        cl += RunWFA2Cpp(sequences, alnP, rounds, loglevel);
        PBLOG_INFO << "WFA2 C++ time : "
                   << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() /
                                                                 rounds / sequences.size());
    }
    if (options[OptionNames::KSW2]) {
        Utility::Stopwatch timer;
        cl += RunKSW2(sequences, alnP, rounds, loglevel);
        PBLOG_INFO << "KSW2 time     : "
                   << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() /
                                                                 rounds / sequences.size());
    }

    PBLOG_DEBUG << cl;

    return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, CreateCLI(), &RunnerSubroutine);
}
