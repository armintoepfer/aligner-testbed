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

using namespace PacBio;

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

int RunnerSubroutine(const CLI_v2::Results& options)
{
    const std::vector<std::string> files = options.PositionalArguments();
    const auto sequences = LoadTargetQueryFromFile(files.front());

    const int32_t seqPairs = sequences.size();
    PBLOG_INFO << "Number of sequence pairs : " << seqPairs;

    Pancake::AlignmentParameters panParameters;

    Utility::Stopwatch timer;
    for (int32_t i = 0; i < 10; ++i) {
        for (const auto& [qry, target] : sequences) {
            wfa::WFAlignerGapAffine2Pieces aligner(
                panParameters.mismatchPenalty, panParameters.gapOpen1, panParameters.gapExtend1,
                panParameters.gapOpen2, panParameters.gapExtend2, wfa::WFAligner::Alignment,
                wfa::WFAligner::MemoryHigh);
            aligner.alignEnd2End(qry.c_str(), qry.size(), target.c_str(), target.size());
        }
    }
    timer.Freeze();
    PBLOG_INFO << "WFA2 time "
               << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() / 10.0 /
                                                             seqPairs);

    timer.Reset();
    auto ksw2Aligner = Pancake::CreateAlignerKSW2(panParameters);
    for (int32_t i = 0; i < 10; ++i) {
        for (const auto& [qry, target] : sequences) {
            ksw2Aligner->Global(qry, target);
        }
    }
    timer.Freeze();
    PBLOG_INFO << "KSW2 time "
               << Utility::Stopwatch::PrettyPrintNanoseconds(timer.ElapsedNanoseconds() / 10.0 /
                                                             seqPairs);

    return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, CreateCLI(), &RunnerSubroutine);
}
