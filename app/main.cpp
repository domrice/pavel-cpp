#include <CLI/CLI.hpp>

int main(int argc, char** argv)
{
    CLI::App app{"PAVEL: decorative patterns with packed volumetric elements"};
    argv = app.ensure_utf8(argv);
    CLI11_PARSE(app, argc, argv);
}
