#ifndef _VERBOSITY_H_
#define _VERBOSITY_H_

namespace chemfem{

  /**
   * Whether the library reports what it does on the console, e.g. the name of every VTK
   * file it writes. On by default. A loop that writes one file per time step or per
   * refinement level should switch it off, its output is more useful than a file name.
   */
  inline bool& VerboseFlag()
  {
    static bool verbose = true;
    return verbose;
  }

  inline void SetVerbose(bool verbose)
  {
    VerboseFlag() = verbose;
  }

  inline bool Verbose()
  {
    return VerboseFlag();
  }

};

#endif
