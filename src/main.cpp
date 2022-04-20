#include <kalloc.h>
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

using namespace PacBio;

namespace OptionNames {
// clang-format off
const CLI_v2::Option WFA2 {
R"({
    "names" : ["wfa2"],
    "description" : "Use WFA2-lib",
    "type" : "bool",
    "default" : true
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
    i.AddOptionGroup("Algorithms", {OptionNames::KSW2, OptionNames::WFA2, OptionNames::MiniWFA});

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

void RunMiniWFA(const std::vector<std::pair<std::string, std::string>>& sequences,
                Pancake::AlignmentParameters& alnP, const int32_t rounds)
{
    Utility::Stopwatch timer;
    mwf_opt_t opt;
    mwf_opt_init(&opt);
    opt.x = alnP.mismatchPenalty;
    opt.o1 = alnP.gapOpen1;
    opt.e1 = alnP.gapExtend1;
    opt.o2 = alnP.gapOpen2;
    opt.e2 = alnP.gapExtend2;
    opt.flag |= MWF_F_CIGAR;
    void* km = 0;

    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            // km = km_init();
            mwf_rst_t rst;
            mwf_wfa(km, &opt, qry.size(), qry.c_str(), target.size(), target.c_str(), &rst);
            // mwf_assert_cigar(&opt, rst.n_cigar, rst.cigar, qry.size(), target.size(), rst.s);
            std::string cigar;
            // for (int32_t i = 0; i < rst.n_cigar; ++i) {
            //     cigar += std::to_string(rst.cigar[i] >> 4) + "MIDNSHP=XBid"[rst.cigar[i] & 0xf];
            // }
            // PBLOG_INFO << rst.s << ' ' << Data::Cigar::FromStdString(cigar).ToStdString();
            kfree(km, rst.cigar);
            // km_destroy(km);
        }
    }

    timer.Freeze();
    PBLOG_INFO << "MWFA time "
               << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() / rounds /
                                                             sequences.size());
}

void RunWFA2(const std::vector<std::pair<std::string, std::string>>& sequences,
             Pancake::AlignmentParameters& alnP, const int32_t rounds)
{
    Utility::Stopwatch timer;
    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            wfa::WFAlignerGapAffine2Pieces aligner(
                alnP.mismatchPenalty, alnP.gapOpen1, alnP.gapExtend1, alnP.gapOpen2,
                alnP.gapExtend2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
            aligner.alignEnd2End(qry.c_str(), qry.size(), target.c_str(), target.size());
            // PBLOG_INFO << aligner.getAlignmentScore() << ' ' << aligner.getAlignmentCigar();
        }
    }
    timer.Freeze();
    PBLOG_INFO << "WFA2 time "
               << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() / rounds /
                                                             sequences.size());
}

void RunKSW2(const std::vector<std::pair<std::string, std::string>>& sequences,
             Pancake::AlignmentParameters& alnP, const int32_t rounds)
{
    Utility::Stopwatch timer;
    auto ksw2Aligner = Pancake::CreateAlignerKSW2(alnP);
    for (int32_t i = 0; i < rounds; ++i) {
        for (const auto& [qry, target] : sequences) {
            const auto res = ksw2Aligner->Global(target, qry);
            // PBLOG_INFO << res.score << ' ' << res.cigar.ToStdString();
        }
    }
    timer.Freeze();
    PBLOG_INFO << "KSW2 time "
               << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() / rounds /
                                                             sequences.size());
}

int RunnerSubroutine(const CLI_v2::Results& options)
{
    const std::vector<std::string> files = options.PositionalArguments();
    const auto sequences = LoadTargetQueryFromFile(files.front());

    const int32_t seqPairs = sequences.size();
    PBLOG_INFO << "Number of sequence pairs : " << seqPairs;

    Pancake::AlignmentParameters alnP;
    const int32_t rounds = options[OptionNames::Rounds];

    if (options[OptionNames::MiniWFA]) {
        RunMiniWFA(sequences, alnP, rounds);
    }
    if (options[OptionNames::WFA2]) {
        RunWFA2(sequences, alnP, rounds);
    }
    if (options[OptionNames::KSW2]) {
        RunKSW2(sequences, alnP, rounds);
    }

    return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, CreateCLI(), &RunnerSubroutine);
}
