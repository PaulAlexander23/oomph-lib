#ifndef PARSE_ARGUMENTS_HEADER
#define PARSE_ARGUMENTS_HEADER

#include "generic.h"

class Arguments
{
public:
  Arguments(int argc, char** argv)
  {
    oomph::CommandLineArgs::setup(argc, argv);
    oomph::CommandLineArgs::specify_command_line_flag(
      "--debug", "Optional: Debug Jacobian. Default = false");
    // Parse and assign command line arguments
    bool has_unrecognised_arg = false;
    oomph::CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);

    debug = oomph::CommandLineArgs::command_line_flag_has_been_set("--debug");
  }

private:
  bool debug;

public:
  bool has_debug()
  {
    return debug;
  }
};


#endif
