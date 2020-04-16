#ifndef AUTOTOOLS_HH
#define AUTOTOOLS_HH

#include <filesystem>

namespace mush
{
namespace fs = std::filesystem;
namespace autotools
{
#ifndef ABS_SRCDIR
#define ABS_SRCDIR "BAD_ABS_SRCDIR_FLAG"
#endif

#ifndef ABS_TOP_BUILDDIR
#define ABS_TOP_BUILDDIR "BAD_ABS_TOP_BUILDDIR_FLAG"
#endif

/// Absolute path to the source directory of the repository, i.e. git root directory
static const fs::path abs_srcdir = ABS_SRCDIR;

/// Absolute path to the build directory, i.e. wherever you ran the configure script
static const fs::path abs_top_builddir = ABS_TOP_BUILDDIR;

/// Absolute path to the input files for tests
static const fs::path input_filesdir = abs_srcdir / "tests/input_files";
/// Absolute path to whatever files get generated during tests
static const fs::path output_filesdir = abs_srcdir / "tests/output_files";

} // namespace autotools
} // namespace mush

#endif
